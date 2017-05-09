#define NMAG 6
#define SIZE 10
#define PART 50

double uniformRandom(void);
double zeroDistTab(void);
int getParameters(void);

char OutFileName[1024];
int ctmax;
int N = 2*SIZE*PART+1;
// arrays of positions
double posx[(2*SIZE*PART+1)*(2*SIZE*PART+1)];
double posy[(2*SIZE*PART+1)*(2*SIZE*PART+1)];
double part_dist[(2*SIZE*PART+1)*(2*SIZE*PART+1)];
//simulation variables
double ek;
int i,j,k;   
double t;   
int run;
double dt; 
double h; 
double g ; 
double mu ; //poczatkowe
double magnetSize;
double abortVel;
double start_amp;
double stop_amp;
double step_amp; 
double start_mu;
double stop_mu;
double step_mu; 
int minSteps;
int maxSteps;
int ct;
int flag;

// Heksagon
 double posxM[NMAG] = {1.732050808,0.0,-1.732050808,1.732050808,0.0,-1.732050808   };
 double posyM[NMAG] = {1.0,2.0,1.0,-1.0,-2.0,-1.0}; 
double min;
// position with minimal distance to magnet
int pos_min;

// distance between magnets and array points
double dist[NMAG];

//velocities
double vpredictx,vpredicty,vx,vy;
//accelerations
double ax,ay,a1x,a1y,a2x,a2y;

double absz;
double rx,ry;
double uk,u0,u;
double dtx,dty;
double last_pos,new_pos;
double u1,u2,u11,u22;	

int white_noise;

double amp;

char buff[80];
double minimalVel = 0.0000001;
