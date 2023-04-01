



void get_rgb(double xyz[3],
              double ulo, double uhi,
              double vlo, double vhi,
              double wlo, double whi,
              int type, double rgb[3]){ // needs vlo,vhi,ulo,uhi for mapping

                double x = 0;
                double y = 0;

      double theta = atan2(xyz[1], xyz[0]);
      double phi = acos(xyz[2]);


      double h_theta = atan2(xyz[2], xyz[0]);
      double h_phi = asinh(xyz[1]);


      if(type == ELLIPSOID){
        x = textureWidth*(theta - ulo)/(uhi - ulo);
        y = textureHeight*(phi - vlo)/(vhi - vlo);
        get_xwd_map_color(idA, x, y, rgb) ; // returns -1 on error, 1 if ok

      }else if(type == HYPER){
        x = textureWidth*(h_theta - ulo)/(uhi - ulo);
        y = textureHeight*(h_phi - wlo)/(whi - wlo);
        get_xwd_map_color(idA, x, y, rgb) ; // returns -1 on error, 1 if ok
      }else if(type == PLANE){
        x = textureWidth*(xyz[2] - ulo)/(uhi - ulo);
        y = textureHeight*(xyz[0] - vlo)/(vhi - vlo);
        get_xwd_map_color(idB, x, y, rgb) ; // returns -1 on error, 1 if ok
      }



}
