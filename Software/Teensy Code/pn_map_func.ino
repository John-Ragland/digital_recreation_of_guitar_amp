void pn_map_func(mtx_type* pn, mtx_type* pn_map){
   // define top and bottom values of pn for given gmap
   int top = 25;
   int bottom = -25;
   // map pn to index of gmap
  
   pn_map[0] = (pn[0] - bottom)*(100/(top-bottom));
   pn_map[1] = (pn[1] - 236)*25;

   if (pn_map[0] <= 1){
    pn_map[0] = 2;
    Serial.println("Table Saturated");
   }
   if (pn_map[0] >= 100){
    pn_map[0] = 100;
    Serial.println("Table Saturated");
   }
   if (pn_map[1] >= 100){
    pn_map[1] = 100;
    Serial.println("Table Saturated");
   }
   if (pn_map[1] >= 100){
    pn_map[1] = 100;
    Serial.println("Table Saturated");
   }                      
}
