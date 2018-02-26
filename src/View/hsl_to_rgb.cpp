#include <math.h>

#include "hsl_to_rgb.hpp"

void hsl_to_rgb( float* rgb,float* hsl){
  float c = (1 - fabs(2*hsl[2] -1))*hsl[1];
  float hp = hsl[0]/60;
  float x= c*(1 - fabs(fmod(hp,2) -1));
  if(hsl[0] != hsl[0] || hp < 0){
    rgb[0]=0;rgb[1]=0;rgb[2]=0;
  }
  else if(hp <= 1){
    rgb[0]=c;rgb[1]=x;rgb[2]=0;
  }
  else if(hp <= 2){
    rgb[0]=x;rgb[1]=c;rgb[2]=0;
  }
  else if(hp <= 3){
    rgb[0]=0;rgb[1]=c;rgb[2]=x;
  }
  else if(hp <= 4){
    rgb[0]=0;rgb[1]=x;rgb[2]=c;
  }
  else if(hp <= 5){
    rgb[0]=x;rgb[1]=0;rgb[2]=c;
  }
  else if(hp <= 6){
    rgb[0]=c;rgb[1]=0;rgb[2]=x;
  }
  float m = hsl[2] - 0.5*c;
  rgb[0] += m; rgb[1] += m; rgb[2] += m;
}

