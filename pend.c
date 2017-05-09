/*

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "pend.h"



int main(int argc, char *argv[]){
	

	if (argc == 1) {
		fprintf(stderr, "\nERROR! No parameters given!\nYou should type: ./pend params\n\n");
		exit(1);
	}

    strcpy (buff, argv[1]);
    getParameters();
    srand (time(NULL));
	zeroDistTab();
	if(!white_noise){
		 start_amp = 0;
		 stop_amp=start_amp;		 
	}
	
	for(mu=start_mu;mu<=stop_mu;mu+=step_mu){
	 for(amp=start_amp;amp<=stop_amp;amp+=step_amp){	
		
	  sprintf(OutFileName, "data-mu_%.3lf-amp_%.1lf",mu,amp);
	  FILE* fp = fopen(OutFileName,"w");
         
         
       for(i=0;i<N;i++){
	    
	    for(j=0;j<N;j++){    
		   
		  posx[i*N+j] = (double)(-SIZE) + (double)1/PART*j;
		  posy[i*N+j] = (double)(-SIZE) + (double)1/PART*i;
			
		     vx =0.0; vy = 0.0; a1x=0.0; a1y=0.0; ax =0.0; ay = 0.0; ct = 0; run =1;
		     
			 while(run == 1){			
				
				last_pos = sqrt(posy[i*N+j]*posy[i*N+j]+posx[i*N+j]*posx[i*N+j]);
				
				dtx = vx*dt+ (4.0*ax - a1x)* dt*dt / 6.0;
				dty =  vy*dt+ (4.0*ay - a1y)* dt*dt / 6.0;  
							
				posy[i*N+j] += dty;
				posx[i*N+j] += dtx;
				
				new_pos = sqrt(posy[i*N+j]*posy[i*N+j]+posx[i*N+j]*posx[i*N+j]);
				
				vpredictx = vx + (3.0*ax - a1x)*dt/2.0;
				vpredicty = vy + (3.0*ay - a1y)*dt/2.0;
								
				a2x = 0.0; a2y = 0.0; uk =0.0; 
				
				   for(k=0;k<NMAG;k++){		 

						if(white_noise){
							u1=uniformRandom();
							u2=uniformRandom();
							//6,555541564 is max value of u11 and u22 because √(−2×(ln(1÷2147483647)))×cos(2×π×0) 
							//2147483647 is RAND_MAX
							u11 = amp*sqrt(-2.0*log(u1))*cos(2.0*3.1415*u2)/6.555541564;  // should be /6,555541564;
							u22 = amp*sqrt(-2.0*log(u1))*sin(2.0*3.1415*u2)/6.555541564;
					    }
					    else {
							u11 = 0; u22=0;
						}
					    absz =  sqrt((posxM[k] - posx[i*N+j])*(posxM[k] - posx[i*N+j]) + (posyM[k] - posy[i*N+j])*(posyM[k] - posy[i*N+j]) );
					
					    a2x += u11+(posxM[k] - posx[i*N+j]) / sqrt((h*h + absz*absz)*(h*h + absz*absz)*(h*h + absz*absz));
					    a2y += u22+(posyM[k] - posy[i*N+j]) / sqrt((h*h + absz*absz)*(h*h + absz*absz)*(h*h + absz*absz)); 			  
						 						 
					    dist[k] =sqrt((posx[i*N+j]-posxM[k])*(posx[i*N+j]-posxM[k]) + ( posy[i*N+j] - posyM[k] )*( posy[i*N+j] - posyM[k] ) );							
					   
					    uk +=  (-1.0)/(sqrt((h*h + absz*absz)*(h*h + absz*absz))); // potentail energy from magnets
					   
					   if(k==0){
						  min = dist[0];
						  pos_min = 1;
					   }
					   if(dist[k] < min){
							min = dist[k];
							pos_min = k+1;										
					   }
					   if(((min < magnetSize)&&(sqrt(vx*vx+vy*vy)<abortVel))) // blisko i mala predkosc
							run = 0;
					   if((sqrt(vx*vx+vy*vy)<minimalVel)&&(ct > minSteps)) // mala predkosc, wykonano duzo krokow, ale dalej niz magnet size
							run = 0;
					   if(ct > maxSteps) // jezeli daleko i duza predkosc przez amplitude white noise
							run = 0;						 
				 }			
	
				 a2x += -g*posx[i*N+j] - mu*vpredictx;
				 a2y += -g*posy[i*N+j] - mu*vpredicty;
				 vx += (2.0*a2x + 5.0*ax - a1x)*dt/6.0;
				 a1x = ax;
				 ax = a2x,
				 vy += (2.0*a2y + 5.0*ay - a1y)*dt/6.0;
				 a1y = ay;
				 ay = a2y;		 			
			     u = g*0.5 *(posx[i*N+j]*posx[i*N+j] + posy[i*N+j]*posy[i*N+j]) + uk;			    
				 ct++; 	
				 
			} //tutaj sie konczy while 					
		    fprintf(fp,"%d ",pos_min);		  
		}	
		 fprintf(fp,"\n");
			
	     if(i%20==0){
			 printf("Row: %d from %d vx: %lf  vy: %lf\n",i,N,vx,vy);
			 printf("Minimal distance: %lf \n",min);
		 }
	
	  }
	  
	  if(white_noise) printf("\nActual amp: %lf \n",amp);
	  fclose(fp);
	 } //loop amp
	} // loop mu
	printf("Data stored in %s\n",OutFileName);
	return 1;
}

int getParameters(void) {
	
FILE *FileIn = fopen (buff, "r"); 
char Line[1024];

 // FILE *FileIn = fopen("params", "r");
  int tmp;
  while (!feof(FileIn)) {
    tmp = fscanf(FileIn, "%s", Line);
    if (feof(FileIn)) continue;
    if (strcmp(Line, "dt") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(dt));
    } else if (strcmp(Line, "h") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(h));
    } else if (strcmp(Line, "g") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(g));
    } else if (strcmp(Line, "magnetSize") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(magnetSize));
    } else if (strcmp(Line, "abortVel") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(abortVel));
    } else if (strcmp(Line, "start_mu") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(start_mu));
    } else if (strcmp(Line, "stop_mu") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(stop_mu));
    }  else if (strcmp(Line, "step_mu") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(step_mu));
    } else if (strcmp(Line, "start_amp") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(start_amp));
    } else if (strcmp(Line, "stop_amp") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(stop_amp));
    }  else if (strcmp(Line, "step_amp") == 0) {
      tmp = fscanf(FileIn, "%lf\n", &(step_amp));
    } else if (strcmp(Line, "minSteps") == 0) {
      tmp = fscanf(FileIn, "%d\n", &(minSteps));
    } else if (strcmp(Line, "maxSteps") == 0) {
      tmp = fscanf(FileIn, "%d\n", &(maxSteps ));
    } else if (strcmp(Line, "white_noise") == 0) {
      tmp = fscanf(FileIn, "%d\n", &(white_noise ));
    }else {
		printf("Unknown parameter: %s\n", Line);
		exit(1);
	  }
  }

  fclose(FileIn);

  
  
  printf("start_mu:\t%.7f\n", start_mu);
  printf("stop_mu:\t%f\n", stop_mu);
  printf("step_mu:\t%f\n", step_mu);
  printf("dt:\t\t%lf\n", dt);  
  printf("start_amp:\t%.7f\n", start_amp);
  printf("stop_amp:\t%f\n", stop_amp);  
  printf("step_amp:\t%f\n", step_amp);
  printf("abortVel:\t%.7f\n", abortVel);
  printf("magnetSize:\t%.7f\n", magnetSize);
  printf("g:\t\t%.7f\n", g);
  printf("h:\t\t%.7f\n", h);
  printf("minSteps:\t%d\n", minSteps);
  printf("maxSteps:\t%d\n", maxSteps);
  printf("white noise:\t%d\n\n", white_noise);

  return 1;
}
double zeroDistTab()
{	
	for(k=0;k<NMAG;k++) 
		dist[k] = 0.0;
	return 1;
}

double uniformRandom()
{	
  return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

