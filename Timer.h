#ifndef CLASS_TIMER
#define CLASS_TIMER

#include <stdio.h>
#include <sys/time.h>
#include <time.h>

using namespace std;

class Timer {
	public:	
		
		typedef unsigned int uint;
		typedef double real;
		bool flag;

		uint c1, c2, cin, cin_tot, cout, cout_tot, cp0, cp;
		real s1, s2, sin, sin_tot, sout, sout_tot, sp0, sp;
	
		Timer(){
			cin=cin_tot=cout=cout_tot=sp0=0; 
			sin=sin_tot=sout=sout_tot=cp0=0;
			flag=false;
			set();
		}
	
		void set(){
			s1= time();
			c1= clock();
			if (flag) {
				sout=s1-s2;
				sout_tot+=sout;
				cout=c1-c2;
				cout_tot+=cout;
			}
			cp=sp=0;
		}
		
		void pause() {
			sp0 = time();
			cp0 = clock();
		}
		
		void resume() {
			if (sp0) { 
				sp+= time() - sp0;
				cp+= clock() - cp0;
				sp0 = cp0 = 0;
			}
		}
	
		double get(){
		
			s2 =  time();
		
			sin=s2-s1-sp;
			sin_tot+=sin;
	
			c2= clock();
			cin=c2-c1-cp;
			cin_tot+=cin;
	
			flag=true;
			
			return sin;
		}
	
		double time() {
			gettimeofday(&tv, 0);
			return tv.tv_sec + tv.tv_usec/1000000.0;
		}
	
		void print(){
			printf("Timer::print()\n");
			printf(" ├ sin  = %g\n",  sin);
			//printf(" ├ sout = %g\n",  sout);
			printf(" └ sin_tot  = %g\n",  sin_tot);
			//printf(" └ sout_tot = %g\n",  sout_tot);
		}
		private:
			timeval tv;

};


#endif
