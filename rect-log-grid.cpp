#include <unistd.h>
#include <dirent.h>
#include <iostream>
#include <sys/stat.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <deque>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <atomic>
#include <thread>

#include "Vector.h"
#include "Timer.h"

 #define SQR(x) ((x)*(x))
using namespace std;
typedef unsigned long int luint;


int e0 = 0, e1 = 4;
uint nn = 16; //32, 16, 8, 4, 2
float t = 0, dt = 1;
uint N = 4;

char dn[] = "/fai/yurin/projects/rect-log-grid/plots-16-004-2";
char fn[256];

float st = .5*pow(2, e0)/nn;

uint nl = e1 - e0 + 1;

uint nn2 = nn*nn;
uint nn3 = nn2*nn;
uint nn3_04 = nn3*04;
uint nn3_16 = nn3*16;
uint nn3_64 = nn3*64;


luint n = luint(nl)*64*nn*nn*nn;
vector<int> vids;
vector<luint> cids;

std::atomic<luint> nops{0};
//luint ops = 0;

luint get_id(luint l, luint i, luint j, luint k, luint ii, luint jj, luint kk) {
	return (((((l*4 + k)*4 + j)*4 + i)*nn + kk)*nn + jj)*nn + ii;
}

struct Cell {
	float s;
	Vector<float> R;
	luint cid;
	int l, i,j,k, ii,jj,kk;
};
vector<Cell> c;

struct Particle {
	Vector<float> R;
	float s; // for plotting purpose only
};
vector<Particle> p;

struct ID {
	int l;
	int i, j, k;
	int ii, jj, kk;
	float s;
	void print(){
		printf("l=%i\n",l);
		printf("i,j,k=%i,%i,%i\n",i,j,k);
		printf("ii,jj,kk=%i,%i,%i\n",ii,jj,kk);
	}
};

Cell get_cell(luint cid) {

	Cell c; uint r;

	c.cid = cid;

	c.l = cid / nn3_64;
	r = cid % nn3_64;

	c.k = r / nn3_16;
	r = r % nn3_16;

	c.j = r / nn3_04;
	r = r % nn3_04;

	c.i = r / nn3;
	r = r % nn3;

	c.kk = r / nn2;
	r  = r % nn2;

	c.jj = r / nn;
	r  = r % nn;

	c.ii = r;

	double s = pow(2, e0+c.l); // WARNING to be optimized later on
	double ss = s/nn;

	c.R =  s*Vector<float>(c.i-2, c.j-2, c.k-2) +
		   ss*Vector<float>(c.ii+.5, c.jj+.5, c.kk+.5);
	c.s = ss;

	return c;

};

ID get_id(Vector<float> R) {

	ID id;

	float x = fabs(R.x);
	float y = fabs(R.y);
	float z = fabs(R.z);

	float xy_max = x < y ? y : x;
	float xyz_max = z < xy_max ? xy_max : z;

	uint e = *(uint*)&xyz_max >> 23;

	if ((int)e - 127 < e0) e = e0 + 127;

	uint uis = e << 23;


	float s = *(float*)&uis;

	float invs = 1/s;

	float ix = R.x*invs+2;
	float iy = R.y*invs+2;
	float iz = R.z*invs+2;

	id.s = s;

	// level
	id.l = e - 127 - e0;

	// base-grid
	id.i = ix;
	id.j = iy;
	id.k = iz;

	// subgrid
	id.ii = (ix - id.i)*nn;
	id.jj = (iy - id.j)*nn;
	id.kk = (iz - id.k)*nn;

	return id;
}

struct Lock {
	std::atomic<bool> flag{false};
	void lock() {
		while ( flag.exchange(true, std::memory_order_relaxed) );
		std::atomic_thread_fence(std::memory_order_acquire);
	}

	void unlock() {
		std::atomic_thread_fence(std::memory_order_release);
		flag.store(false, std::memory_order_relaxed);
	}
};

std::deque< Lock > locks;
	

void plot_grid(const char* fn) {

	FILE *gp = popen("gnuplot --persist", "w");
	//fprintf(gp, "set term wxt size 800, 800 font 'Ubuntu, 7'\n");


	double W = 2*1920, H = 2*1080;
	double lm=20, rm=10, tm=10, bm=18;

	fprintf(gp, "set term pngcairo size %g, %g font 'Ubuntu, 7' enhanced\n", W, H);
	fprintf(gp, "set output '%s'\n", fn);


	fprintf(gp, "set tmargin at screen %g\n", 1-tm/H);
	fprintf(gp, "set bmargin at screen %g\n", bm/H);
	fprintf(gp, "set lmargin at screen %g\n", lm/W);
	fprintf(gp, "set rmargin at screen %g\n", 1-rm/W);

	fprintf(gp, "set colors classic\n");
	fprintf(gp, "set xtics 1\n");
	fprintf(gp, "set ytics 1\n");
	fprintf(gp, "set grid\n");

	float min =  1;//*min_element(vids.begin(), vids.end());
	float max =  p.size();//*max_element(vids.begin(), vids.end());

	fprintf(gp, "set cbrange [%g:%g]\n", min, max);
	//fprintf(gp, "set logscale cb\n", min, max);
	//fprintf(gp, "set palette defined ( 0 'black', 1 'blue', 2 'red', 3 'orange')\n");
	fprintf(gp, "set palette defined ( 0 'dark-green',1 'green', 2 'blue', 3 'dark-blue', 4 'red', 5 'dark-red', 6 'orange', 7 'black' )\n");

	fprintf(gp, "set label 1 \"t = %g\\ndt = %g\\nncells = %i (active)\\nncells = %lu (total)\\nnops =  %lu\\nnn = %i\\nnparticles = %i\\nrange = %g - %g\\nnthreads = %i\" at screen %g, %g font 'Arial, 18'\n", t, dt, cids.size(), vids.size(), (luint)nops, nn, p.size(), pow(2,e0), pow(2,e1), thread::hardware_concurrency(), 2*lm/W, 1-3*tm/H);


	uint id = 2;
	for (int l=0; l<e1-e0+1; l++)
		for (uint j=0; j<4; j++)
		for (uint i=0; i<4; i++)
		if ( (j!=1 && j!=2) ||
			  (i!=1 && i!=2) || l==0 ) {

			double s = pow(2, e0+l);
			double x = (i - 1.5)*s;
			double y = (j - 1.5)*s;

			fprintf(gp, "set object %i rectangle center %g, %g size %g,%g fs empty border 0 front\n", ++id, x, y, s,s);

			for (uint jj=0; jj<nn; jj++)
			for (uint ii=0; ii<nn; ii++) {

				double ss = s/nn;
				double xx = x -.5*s + (ii + 0.5)*ss;
				double yy = y -.5*s + (jj + 0.5)*ss;

				float v = vids[get_id(l,i,j,2,ii,jj,0)];

				fprintf(gp, "set object %i rectangle center %g, %g size %g,%g fc palette cb %g fs solid %g border rgb '#%x' lw 0.01\n", ++id, xx, yy, ss,ss, v, v ? 0.5 : 0, v ? 0x444444: 0xBBBBBB);
			}
		}

	for (auto& cid:cids) {
		Cell c = get_cell(cid);
		if (c.R.z-.5*c.s == 0)
		fprintf(gp, "set object %i rectangle center %g, %g size %g,%g  fc palette cb %g fs solid %g border rgb '#%x' lw 0.01\n", 
				  ++id, c.R.x, c.R.y, c.s,c.s, (float)vids[cid], vids[cid] ? 0.7 : 0, vids[cid] ? 0x111111: 0xBBBBBB );
	}

	
	for (int j=0; j<p.size(); j++) { float f = 2;
		if (t<f*p[j].s)
		fprintf(gp, "set object %i circle at %g, %g size %g fs empty border lc palette cb %i lw %g\n", 
				  ++id, p[j].R.x, p[j].R.y, t, j+1, t<p[j].s ? 1-t/(f*p[j].s) : 0.); }
	
	

	double s = 16;
	double asp = (W-lm-rm)/(H-tm-bm);

	fprintf(gp, "plot [-%g:+%g][-%g:+%g]", asp*s, asp*s, s,s);
	fprintf(gp, " '-' u 1:2:3 w p pt 7 ps 0.7 palette t ''");
	//fprintf(gp, ",'-' u 1:2:3:4 w circles lt 1 palette lw var t ''");
	fprintf(gp,"\n");
	
	for (int j=0; j<p.size(); j++)
	fprintf(gp, "%g %g %i\n", p[j].R.x, p[j].R.y, j+1);
	fprintf(gp, "e\n");

	/*
	for (int j=0; j<p.size(); j++) {

		//for (auto& cid:cids) if (vids[cid]==j+1) {
		//if 
		if (t < .5*p[j].s )
		fprintf(gp, "%g %g %g %i 0.5\n", p[j].R.x, p[j].R.y, t, j+1 );  
		//break; }

	}
	fprintf(gp, "e\n");
	/**/
	
	fclose(gp);
}


void init_rad_transfer() {

	for (int j=0; j<p.size(); j++) {

		ID id = get_id(p[j].R);
		luint cid = get_id(id.l, id.i,id.j,id.k, id.ii,id.jj,id.kk);
		cids.push_back(cid);
		vids[cid]=j+1;

		//printf("%i  %i %i %i  %2i %2i %2i\n", id.l, id.i,id.j,id.k, id.ii, id.jj, id.kk);
	}

	//printf("cids.size() = %i\n\n", st, cids.size());
}

void iterate_rad_transfer() {
	
	vector<thread> th(thread::hardware_concurrency());
	//th.resize(1);
	if (cids.size()<th.size()) th.resize(cids.size());


	vector< vector<luint> > cids_new(th.size());

	vector<Timer> tr(th.size());
	//vector<luint> nops(th.size());

	uint nth = ceil(cids.size()/(double)th.size());

	for (uint k=0; k<th.size(); k++) th[k] =
			thread([&](uint k, uint i0, uint i1) {


		vector<Cell> cs; cs.reserve(4*1024);
		cids_new[k].reserve(4*1024);

		for (int i=i0; i<i1; i++) {

			Cell c = get_cell(cids[i]);

			//printf("s = %g, t = %g   %i\n", c.s/st, t/st,  uint(t/st)%uint(c.s/st) );
			//float _cs = st<c.s ? .5*c.s : c.s; // we update each cell at
			if (uint(t/st)%uint(.5*c.s/st)!=0) {
				cids_new[k].push_back(cids[i]);
				continue;
			}

			cs.clear();
			//printf("%7lu: %2i %2i %2i  % 7.3f % 7.3f % 7.3f  %g\n", cids[i], c.ii, c.jj, c.kk, c.R.x,  c.R.y, c.R.z, c.s );

			for (int kk=c.kk-1; kk<=c.kk+1; kk++)
			for (int jj=c.jj-1; jj<=c.jj+1; jj++)
			for (int ii=c.ii-1; ii<=c.ii+1; ii++)

			if ( ii<0 || nn<=ii ||
				  jj<0 || nn<=jj ||
				  kk<0 || nn<=kk ) {

				Vector<> R = c.R + c.s*Vector<float>(ii-c.ii, jj-c.jj, kk-c.kk);

				if (e1<=e0+c.l)
				if ( 2*c.s*nn<fabs(R.x) ||
					  2*c.s*nn<fabs(R.y) ||
					  2*c.s*nn<fabs(R.z) ) continue;

				std::unordered_set<luint> uset;

				for (int kkk=0; kkk<2; kkk++)
				for (int jjj=0; jjj<2; jjj++)
				for (int iii=0; iii<2; iii++) {

					Vector<float> RR = R + .5*c.s*Vector<float>(iii-.5, jjj-.5, kkk-.5);
					ID id = get_id(RR);
					uset.insert( get_id(id.l, id.i,id.j,id.k, id.ii,id.jj,id.kk) );
					//cs.push_back(get_cell(get_id(id.l, id.i,id.j,id.k, id.ii,id.jj,id.kk)));

				}

				for (auto& us:uset)
					cs.push_back(get_cell(us));

			} else
				cs.push_back(get_cell(get_id(c.l, c.i,c.j,c.k, ii,jj,kk )));


			bool freeCellsF = false;

			for (auto& c:cs) {

				locks[c.cid].lock();

				if (vids[c.cid]==0) {

					if ((c.R-p[vids[cids[i]]-1].R).sqr()<t*t) {
						vids[c.cid] = vids[cids[i]];
						cids_new[k].push_back(c.cid);
						nops++;
					} else  freeCellsF = true;
					
				} else if (vids[c.cid]!=vids[cids[i]]) {

					float r2_0 = (c.R-p[vids[cids[i]]-1].R).sqr();
					float r2_1 = (c.R-p[vids[c.cid]-1].R).sqr();

					if (r2_0<r2_1) vids[c.cid] = vids[cids[i]];
					nops++;
				}

				locks[c.cid].unlock();
			}

			if (freeCellsF) cids_new[k].push_back(cids[i]);
		}

	}, k, (k+0)*nth, (k+1)*nth<cids.size() ? (k+1)*nth : cids.size() );

	for (auto& thi : th) thi.join();

	cids.clear();
	for (int k=0; k<th.size(); k++)
		for (int i=0; i<cids_new[k].size(); i++)
			cids.push_back(cids_new[k][i]);

}

float get_new_dt() {

	float ss_min = 1e30;
	for (int i=0; i<cids.size(); i++) {
		Cell c = get_cell(cids[i]);
		if (.5*c.s<ss_min) ss_min = .5*c.s;
	}

	if (uint(t/st)%uint(ss_min/st)==0) return ss_min;
	else return dt;
}


void mark_voronoi_direct() {

	vector<thread> th(thread::hardware_concurrency());
	//th.resize(1);

	luint nth = ceil(vids.size()/(double)th.size());

	for (uint k=0; k<th.size(); k++) th[k] =
		thread([&](uint k, luint cid0, luint cid1) {

		for (luint cid=cid0; cid<cid1; cid++) {

			Cell c = get_cell(cid);

			if ( c.k<1 || 2<c.k ||
				  c.j<1 || 2<c.j ||
				  c.i<1 || 2<c.i || c.l==0 ) {

				float r2_min = 1e43; uint j_min;

				for (int j=0; j<N; j++) {
					float r2 = (c.R-p[j].R).sqr();
					if ( r2 < r2_min) { r2_min = r2; j_min = j; }
				}

				vids[cid] = j_min+1;
			}

		}

	}, k, (k+0)*nth, (k+1)*nth<vids.size() ? (k+1)*nth : vids.size() );

	for (auto& thi : th) thi.join();

}

void mark_voronoi_radial() {

	Timer tr;
	tr.set();
	locks.resize(n);
	printf("locks allocation of %lu elements (or %g MB) took %g secs\n", n,  1*n/(1024*1024.), tr.get());

	init_rad_transfer();

	uint step(0);

	do {

		printf("%3i: t = %6.4f, dt = %6.4f, cids.size() = %-10i  %lu\n", step, t, dt, cids.size(), (luint)nops);
		iterate_rad_transfer();

		//if (i%1==0) { sprintf(fn, "%s/grid.%03i.png", dn, i+1); plot_grid(fn); }
		if (cids.size()==0) break;

		dt = get_new_dt();
		t+=dt;
		step++;

	} while (true);

}


void make_voronoi_radial_marking_plot() {

	Timer tr;
	tr.set();
	locks.resize(n);
	tr.get(); tr.print();

	for (int j=0; j<N; j++) {
		ID id = get_id(p[j].R);
		p[j].s = id.s;
	}

	init_rad_transfer();

	uint step(0);

	do {

		printf("%3i: t = %6.4f, dt = %6.4f, cids.size() = %-10i  %lu\n", step, t, dt, cids.size(), (luint)nops);
		iterate_rad_transfer();

		if (step%1==0) { sprintf(fn, "%s/grid.%03i.png", dn, step+1); plot_grid(fn); }
		if (cids.size()==0) break;

		dt = get_new_dt();
		t+=dt;
		step++;

	} while (true);

}



int main(int argc, char *argv[]) {

	
	Timer tr; tr.set();

	vids.resize(n);

	printf("vids allocation of %lu elements (or %g MB) took %g secs\n", n,  4*n/(1024*1024.), tr.get());
	printf("range: %g - %g\n", pow(2,e0), pow(2,e1));

	srand48(4);

	p.resize(N);
	for (int j=0; j<N; j++) {// WARNING particles has to be in [2^e0, 2^e1] range!
		p[j].R = 32*drand48()*Vector<float>().rnd(); p[j].R.z = 0;
		//32*Vector<float>(1-2*drand48(), 1-2*drand48(), 0*(1-2*drand48()));
	}



	make_voronoi_radial_marking_plot(); exit(0);

	tr.set();
	mark_voronoi_radial();
	printf("radial method took %g secs\n", tr.get());


	tr.set();
	mark_voronoi_direct();
	printf("direct method took %g secs\n", tr.get());


	printf("th.size(): %lu\n", thread::hardware_concurrency());
	printf("vids.size(): %lu\n", vids.size());
	printf("p.size(): %u\n", p.size());


	return 0;
	
}
