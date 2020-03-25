#include <MatrixMath.h>
#include "Digital_Amp_vars.h"
//#include "Digital_Amp_funcs.h"
#include <TimerOne.h>
#include <SPI.h>
#include<CircularBuffer.h>

// Filter Clock Pin Number
const int filt_clk = 0;

// Chip Select Pin Numbers
const int CS_ADC = 10;
const int CS_DAC = 9;

// Analog Value (16 bits) and Float
uint16_t analog_value_int = 0;
float analog_value_flt = 0;

uint16_t dac_out;

// Used to Debug Timing
const int test_pin = 2;

// Values of Control Pots
int bass_int, mid_int, treb_int, vol_int, rev_int, vol_discrete_int;
float bass, mid, treb, vol, rev, vol_discrete;


// Input and Output Audio Buffers
CircularBuffer<float,4> aud_in;
CircularBuffer<float,4> aud_out;
CircularBuffer<float,4> toneStack_in;
CircularBuffer<float,4> toneStack_out;

// Calulated output value
float calc_output_flt;
uint16_t calc_output_int;

// Tube simulation IO variables
float tube_pre_in;
float tube_pre_out;
float tube_pre_out2;


double Az[4] = {1,  0, 0, 0};
double Bz[4] = {1, 0, 0, 0};

void setup() {
  Serial.begin(38400);

  // set timing test pin to output
  pinMode(test_pin, OUTPUT);
  
  // Define ADC CS Pin
  pinMode(CS_ADC, OUTPUT);
  digitalWrite(CS_ADC,HIGH);

  // Define DAC CS Pin
  pinMode(CS_DAC, OUTPUT);
  digitalWrite(CS_DAC,HIGH);

  // Begin SPI
  SPI.begin();

  // Set Up Elliptical Filter Sample Frequency (1 MHz)
  pinMode(filt_clk,OUTPUT);
  analogWriteFrequency(filt_clk, 1000000);
  analogWrite(filt_clk,127);

  // Initialize Audio Sample Interrupt to ~ 44.1 kHz
  Timer1.initialize(23);
  Timer1.attachInterrupt(sample_interrupt);
}

void loop() {
  // Read Pot Values
  bass_int = analogRead(0);
  mid_int = analogRead(1);
  treb_int = analogRead(2);
  vol_int = analogRead(3);
  rev_int = analogRead(4);

  bass = (float) bass_int/1024;
  mid = (float) mid_int/1024;
  treb = (float) treb_int/1024;
  vol = (float) vol_int/1024;
  rev = (float) rev_int/1024;
 
  bass = 1 - bass;
  mid = 1 - mid;
  treb = 1 - treb;
  vol = 1 - vol;
  rev = 1 - rev;

  // Compute Transfer Function from Pot Values
  get_tf(bass, mid, treb);

  get_sim_mats();
  

  // Print Pot Values for Debug
  /*
  Serial.print(bass);
  Serial.print("  ");
  Serial.print(mid);
  Serial.print("  ");
  Serial.print(treb);
  Serial.print("  ");
  Serial.print(vol);
  Serial.print("  ");
  Serial.println(rev);
  */
  
  delay(10);

}

void sample_interrupt(void) {
  
  // write previous output to DAC
  DAC_Write(dac_out);
  Serial.print(dac_out);
  Serial.print("  ");
  // read analog value
  analog_value_int = ADC_Read();
  
  // Converts to Input Voltage Values
  analog_value_flt = ((float) analog_value_int - 27218)*0.2/23120;
  aud_in.unshift(analog_value_flt);
  Serial.println(analog_value_int);

 
  // Simulate 2nd Stage of Pre-Amp Circuit
  tube_pre_in = aud_in[0];
  tube_pre_amp();
  
  // quick (inaccurate) scaling
  tube_pre_out = (tube_pre_out + 500000000)/1000000000;
  
  
  // Set Input to Tone Stack
  toneStack_in.unshift(tube_pre_out);
  // Calculate Tone Stack Output
  digitalWrite(test_pin,HIGH);
  toneStack();

  
  digitalWrite(test_pin,LOW);
  // Set Output Audio to Tone Stack Output
  aud_out.unshift(toneStack_out[0]);

  // Scales float output to int DAC output - Is variable but depends on scaling from ADC
  /* commented for quick and dirty scaling
  dac_out = (uint16_t)((aud_out[0]*23120/0.2) + 27218);
  */  

  dac_out = (uint16_t)((aud_out[0]*65536) + 32768);

}

uint16_t ADC_Read(void)
{
  SPISettings SPI_settings(4700000, MSBFIRST, SPI_MODE0);
  SPI.beginTransaction(SPI_settings);
  digitalWrite(CS_ADC, LOW);  // Activate Chip Select
  
  // there is no MOSI for this chip so lets just send bytes of 0s
  byte leading8zeros = SPI.transfer(0x00);
  byte hi = SPI.transfer(0x00);
  byte lo = SPI.transfer(0x00);
  digitalWrite(CS_ADC, HIGH); //deactivate Chipselect

  SPI.endTransaction();
  
  uint16_t adcValue = (hi << 8) + lo; // combining the 2 return Values
  return adcValue;
}

void DAC_Write(uint16_t analog_value)
{
  uint8_t val_low = analog_value & 0xff;
  uint8_t val_high = analog_value >> 8;

  SPISettings SPI_settings(10000000, MSBFIRST, SPI_MODE3);
  SPI.beginTransaction(SPI_settings);
  digitalWrite(CS_DAC, LOW);  // Activate Chip Select

  SPI.transfer(val_high);
  SPI.transfer(val_low);

  digitalWrite(CS_DAC, HIGH);
  SPI.endTransaction();
 
}

void get_tf(float bass, float mid, float treb)
{
  
  float R1 = 250e+3;
  float R2 = 250e+3;
  float R3 = 25e+3;
  float R4 = 100e+3;
  float C1 = 0.250e-9;
  float C2 = 22e-9;
  float C3 = 22e-9;
  float b1, b2, b3, a0, a1, a2, a3;
  float Bz0, Bz1, Bz2, Bz3, Az0, Az1, Az2, Az3;
  b1 = treb*C1*R1 + mid*C3*R3+bass*(C1*R2+C2*R2)+(C1*R3 + C2*R3);

  b2 = treb*(C1*C2*R1*R4+C1*C3*R1*R4)-pow(mid,2)*(C1*C3*pow(R3,2)+C2*C3*pow(R3,2))
    + mid*(C1*C3*R1*R3+C1*C3*pow(R3,2)+C2*C3*pow(R3,2))
    + bass*(C1*C2*R1*R2 + C1*C2*R2*R4+C1*C3*R2*R4)
    + bass*mid*(C1*C3*R2*R3+C2*C3*R2*R3)
    + (C1*C2*R1*R3+C1*C2*R3*R4+C1*C3*R3*R4);

  b3 = bass*mid*(C1*C2*C3*R1*R2*R3+C1*C2*C3*R2*R3*R4)
    -pow(mid,2)*(C1*C2*C3*R1*pow(R3,2)+C1*C2*C3*pow(R3,2)*R4)
    +mid*(C1*C2*C3*R1*pow(R3,2)+C1*C2*C3*pow(R3,2)*R4)
    +treb*C1*C2*C3*R1*R3*R4-treb*mid*C1*C2*C3*R1*R3*R4
    +treb*bass*C1*C2*C3*R1*R2*R4;

  a0 = 1;

  a1 = (C1*R1+C1*R3+C2*R3+C2*R4+C3*R4)+mid*C3*R3+bass*(C1*R2+C2*R2);

  a2 = mid*(C1*C3*R1*R3-C2*C3*R3*R4+C1*C3*pow(R3,2) 
    + C2*C3*pow(R3,2))+bass*mid*(C1*C3*R2*R3+C2*C3*R2*R3) 
    - pow(mid,2)*(C1*C3*pow(R3,2)+C2*C3*pow(R3,2))+bass*(C1*C2*R2*R4 
    + C1*C2*R1*R2+C1*C3*R2*R4+C2*C3*R2*R4)
    +(C1*C2*R1*R4+C1*C3*R1*R4+C1*C2*R3*R4
    +C1*C2*R1*R3+C1*C3*R3*R4+C2*C3*R3*R4);

  a3 = bass*mid*(C1*C2*C3*R1*R2*R3+C1*C2*C3*R2*R3*R4)
    -pow(mid,2)*(C1*C2*C3*R1*pow(R3,2)+C1*C2*C3*pow(R3,2)*R4)
    +mid*(C1*C2*C3*pow(R3,2)*R4+C1*C2*C3*R1*pow(R3,2) 
    -C1*C2*C3*R1*R3*R4) + bass*C1*C2*C3*R1*R2*R4 
    + C1*C2*C3*R1*R3*R4;


  // Creating of Z Domain Filter with Bilinear Transform
  Bz0 = -b1*c - b2*pow(c,2) - b3*pow(c,3);
  Bz1 = -b1*c + b2*pow(c,2) + 3*b3*pow(c,3);
  Bz2 = b1*c + b2*pow(c,2) - 3*b3*pow(c,3);
  Bz3 = b1*c - b2*pow(c,2) + b3*pow(c,3);

  Az0 = -a0 - a1*c - a2*pow(c,2) - a3*pow(c,3);
  Az1 = -3*a0 - a1*c + a2*pow(c,2) + 3*a3*pow(c,3);
  Az2 = -3*a0 + a1*c + a2*pow(c,2) - 3*a3*pow(c,3);
  Az3 = -a0 + a1*c - a2*pow(c,2) + a3*pow(c,3);

  Bz[0] = Bz0;
  Bz[1] = Bz1;
  Bz[2] = Bz2;
  Bz[3] = Bz3;
  Az[0] = Az0;
  Az[1] = Az1;
  Az[2] = Az2;
  Az[3] = Az3;

/*
  Serial.print(Bz0);
  Serial.print(' ');
    Serial.print(Bz1);
  Serial.print(' ');
    Serial.print(Bz2);
  Serial.print(' ');
    Serial.print(Bz3);
  Serial.print(' ');
    Serial.print(Az0);
  Serial.print(' ');
      Serial.print(Az1);
  Serial.print(' ');
      Serial.print(Az2);
  Serial.print(' ');
      Serial.print(Az3);
  Serial.println(' ');
*/
   
}

void toneStack(){
  // Output is returned as global variable
  calc_output_flt = ((toneStack_in[0]*Bz[0] + toneStack_in[1]*Bz[1] + toneStack_in[2]*Bz[2] + toneStack_in[3]*Bz[3]) - (toneStack_out[0]*Az[1] + toneStack_out[1]*Az[2] + toneStack_out[2]*Az[3]))/Az[0];
  toneStack_out.unshift(calc_output_flt);
  
}

void tube_pre_amp(){
  un[1][0] = 51.2*aud_in[0];


// Calculate Pn
  Matrix.MultiplyScalar((mtx_type*) D,(float) Fs,2,2,(mtx_type*)temp22);
  Matrix.Multiply((mtx_type*) temp22,(mtx_type*) H,2,2,2,(mtx_type*)temp22_2);
  Matrix.Multiply((mtx_type*) temp22_2, (mtx_type*) xnm1,2,2,1,(mtx_type*)temp21);



  Matrix.Multiply((mtx_type*) D,(mtx_type*) H,2,2,2,(mtx_type*)temp22);
  Matrix.Multiply((mtx_type*) temp22,(mtx_type*) B,2,2,2,(mtx_type*)temp22_2);
  Matrix.Add((mtx_type*)temp22_2,(mtx_type*)E,2,2,(mtx_type*)temp22);
  Matrix.Multiply((mtx_type*) temp22,(mtx_type*) un,2,2,1,(mtx_type*)temp21_2);
  Matrix.Add((mtx_type*)temp21,(mtx_type*)temp21_2,2,1,(mtx_type*)pn);


  
// Calculate Pn Mapping for g_map
  pn_map_func((mtx_type*) pn, (mtx_type*)pn_map);

  //round index up to next integer
  pn_map_round[0][0] = (int) round(pn_map[0][0]+0.5);
  pn_map_round[1][0] = (int) round(pn_map[1][0]+0.5);

  //Create Index Values for Bilinear Interpolation
  x1 = (pn_map_round[0][0] - 1);
  x2 = (pn_map_round[0][0]);
  y_1 = (pn_map_round[1][0] - 1);
  y2 = (pn_map_round[1][0]);
  x = pn_map[0][0];
  y = pn_map[1][0];

  // Compute Bilinear Interpolation
  Q11 = g_map1[x1][y_1];
  Q12 = g_map1[x1][y2];
  Q21 = g_map1[x2][y_1];
  Q22 = g_map1[x2][y2];

  R1 = ((x2 - x)/(x2 - x1))*Q11 + ((x - x1)/(x2 - x1))*Q21;
  R2 = ((x2 - x)/(x2 - x1))*Q12 + ((x - x1)/(x2 - x1))*Q22;
  in1 = ((y2 - y)/(y2 - y_1))*R1 + ((y - y_1)/(y2 - y_1))*R2;
  
  Q11 = g_map2[x1][y_1];
  Q12 = g_map2[x1][y2];
  Q21 = g_map2[x2][y_1];
  Q22 = g_map2[x2][y2];

  R1 = ((x2 - x)/(x2 - x1))*Q11 + ((x - x1)/(x2 - x1))*Q21;
  R2 = ((x2 - x)/(x2 - x1))*Q12 + ((x - x1)/(x2 - x1))*Q22;
  in2 = ((y2 - y)/(y2 - y_1))*R1 + ((y - y_1)/(y2 - y_1))*R2;

  in[0][0] = in1;
  in[1][0] = in2;
  
// Calculate Xn
  Matrix.MultiplyScalar((mtx_type*)H,(float)Fs,2,2,(mtx_type*)temp22);
  Matrix.Multiply((mtx_type*)temp22,(mtx_type*)xnm1,2,2,1,(mtx_type*)temp21);
  
  Matrix.Multiply((mtx_type*)B,(mtx_type*)un,2,2,1,(mtx_type*)temp21_2);
  Matrix.Multiply((mtx_type*)C,(mtx_type*)in,2,2,1,(mtx_type*)temp21_3);
  Matrix.Add((mtx_type*)temp21_2,(mtx_type*)temp21_3,2,1,(mtx_type*)temp21_4);
  Matrix.Multiply((mtx_type*)H,(mtx_type*)temp21_4,2,2,1,(mtx_type*)temp21_2);
  Matrix.Add((mtx_type*)temp21,(mtx_type*)temp21_2,2,1,(mtx_type*)xn);
  //Matrix.Print((mtx_type*) xn,2,1,"xn");

// Calculate Output
  Matrix.Multiply((mtx_type*)L,(mtx_type*)xn,1,2,1,(mtx_type*)temp21);
  Matrix.Multiply((mtx_type*)M,(mtx_type*)un,1,2,1,(mtx_type*)temp21_2);
  Matrix.Multiply((mtx_type*)N,(mtx_type*)in,1,2,1,(mtx_type*)temp21_3);
  Matrix.Add((mtx_type*)temp21,(mtx_type*)temp21_2,1,1,(mtx_type*)temp21_4);
  Matrix.Add((mtx_type*)temp21_4,(mtx_type*)temp21_3,1,1,(mtx_type*)out);
  //Matrix.Print((mtx_type*)out,1,1,"output");

  tube_pre_out = out[0][0];

}

void get_sim_mats(){
  vol_discrete = round(vol*20);

  vol_discrete_int = (int) vol_discrete;


  
  Serial.print("  ");
  
  if(vol_discrete_int == 0){
      Matrix.Copy((mtx_type*)D_0,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_0,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_0,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_0,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_0,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_0,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_0,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_0,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_0,2,2,(mtx_type*)K);
  }
  if(vol_discrete_int == 1){
      Matrix.Copy((mtx_type*)D0_5,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H0_5,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B0_5,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C0_5,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E0_5,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L0_5,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M0_5,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N0_5,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K0_5,2,2,(mtx_type*)K);
        }
  if(vol_discrete_int == 2){
      Matrix.Copy((mtx_type*)D_01,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_01,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_01,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_01,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_01,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_01,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_01,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_01,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_01,2,2,(mtx_type*)K);
  }
  if(vol_discrete_int == 3){
      Matrix.Copy((mtx_type*)D_15,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_15,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_15,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_15,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_15,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_15,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_15,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_15,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_15,2,2,(mtx_type*)K);
      }

  if(vol_discrete_int == 4){
      Matrix.Copy((mtx_type*)D_02,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_02,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_02,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_02,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_02,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_02,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_02,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_02,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_02,2,2,(mtx_type*)K);
    }

  if(vol_discrete_int == 5){
      Matrix.Copy((mtx_type*)D_25,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_25,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_25,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_25,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_25,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_25,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_25,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_25,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_25,2,2,(mtx_type*)K);
        }

  if(vol_discrete_int == 6){
      Matrix.Copy((mtx_type*)D_03,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_03,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_03,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_03,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_03,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_03,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_03,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_03,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_03,2,2,(mtx_type*)K);
  }
  if(vol_discrete_int == 7){
      Matrix.Copy((mtx_type*)D_35,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_35,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_35,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_35,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_35,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_35,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_35,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_35,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_35,2,2,(mtx_type*)K);
  }
  if(vol_discrete_int == 8){
      Matrix.Copy((mtx_type*)D_04,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_04,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_04,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_04,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_04,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_04,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_04,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_04,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_04,2,2,(mtx_type*)K);
  }
  if(vol_discrete_int == 9){
      Matrix.Copy((mtx_type*)D_45,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_45,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_45,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_45,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_45,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_45,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_45,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_45,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_45,2,2,(mtx_type*)K);
  }

  if(vol_discrete_int == 10){
      Matrix.Copy((mtx_type*)D_05,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_05,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_05,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_05,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_05,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_05,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_05,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_05,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_05,2,2,(mtx_type*)K);
  }

  if(vol_discrete_int == 11){
      Matrix.Copy((mtx_type*)D_55,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_55,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_55,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_55,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_55,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_55,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_55,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_55,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_55,2,2,(mtx_type*)K);
  }

  if(vol_discrete_int == 12){
      Matrix.Copy((mtx_type*)D_06,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_06,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_06,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_06,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_06,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_06,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_06,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_06,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_06,2,2,(mtx_type*)K);
  }
  if(vol_discrete_int == 13){
      Matrix.Copy((mtx_type*)D_65,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_65,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_65,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_65,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_65,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_65,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_65,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_65,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_65,2,2,(mtx_type*)K);
  }
  if(vol_discrete_int == 14){
      Matrix.Copy((mtx_type*)D_07,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_07,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_07,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_07,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_07,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_07,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_07,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_07,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_07,2,2,(mtx_type*)K);
  }
  if(vol_discrete_int == 15){
      Matrix.Copy((mtx_type*)D_75,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_75,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_75,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_75,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_75,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_75,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_75,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_75,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_75,2,2,(mtx_type*)K);
  }
  if(vol_discrete_int == 16){
      Matrix.Copy((mtx_type*)D_08,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_08,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_08,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_08,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_08,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_08,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_08,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_08,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_08,2,2,(mtx_type*)K);
  }
  
  if(vol_discrete_int == 17){
      Matrix.Copy((mtx_type*)D_85,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_85,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_85,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_85,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_85,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_85,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_85,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_85,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_85,2,2,(mtx_type*)K);
  }
  
  if(vol_discrete_int == 18){
      Matrix.Copy((mtx_type*)D_09,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_09,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_09,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_09,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_09,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_09,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_09,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_09,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_09,2,2,(mtx_type*)K);
  }

  if(vol_discrete_int == 19){
      Matrix.Copy((mtx_type*)D_95,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_95,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_95,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_95,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_95,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_95,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_95,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_95,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_95,2,2,(mtx_type*)K);
  }

  if(vol_discrete_int == 20){
      Matrix.Copy((mtx_type*)D_10,2,2,(mtx_type*)D);
      Matrix.Copy((mtx_type*)H_10,2,2,(mtx_type*)H);
      Matrix.Copy((mtx_type*)B_10,2,2,(mtx_type*)B);
      Matrix.Copy((mtx_type*)C_10,2,2,(mtx_type*)C);
      Matrix.Copy((mtx_type*)E_10,1,2,(mtx_type*)E);
      Matrix.Copy((mtx_type*)L_10,1,2,(mtx_type*)L);
      Matrix.Copy((mtx_type*)M_10,1,2,(mtx_type*)M);
      Matrix.Copy((mtx_type*)N_10,1,2,(mtx_type*)N);
      Matrix.Copy((mtx_type*)K_10,2,2,(mtx_type*)K);
  }

}
