/*
1. �� ������� ������� ������ � �������� �������, �� ��������������.
��� �������� � ������ ���� ������ ������.
*/

/*
������������� ������
*/
#include "stdafx.h"
#include <time.h>
#include <fstream>
#include <math.h>
#include <cmath>
//#include <vector>
//#include <exception>

//using namespace System;

#define _CRT_SECURE_NO_WARNINGS

/*
���������� ���������
*/

// ����� ������� ���������� ������, �� ������� ����������� ����� 
short K, K2;

// ����� ������ �� ��� ������
int NP = 6976;

int COLL_COUNT = 0;

// delete this dirty string
#define PI 3.141592653589793238462

// ��������� ������, �������� � load()
double A, A2, dA, L, dL, global_E = 0.0;

// ������ ���������� �������� � ������� �������.
int last;

// ������ "�������"
typedef struct Event_ {
	double t;
	int im, jm;
} Event;

// ������� ������� - ���������� 8192*2 ��������
Event time_queue[16384];

// ������ "�������"
// x,y,z,vx,vy,vz,t - particle coordinats.
// dt - delta time for the next event of this particle.
// kv - velocity koeficient.
typedef struct particle_ {
	double x, y, z, vx, vy, vz, t, dt, kv;
	int x_box, y_box, z_box, ti, box_i, i_copy;
} particle;

// ������ ������, ������ ������� N*2 + ����������
// � ������� ������� � ����� ������ ������� ������.
particle particles[16384];

// ������. ����� ������� ������� �� ��������� ������,
// ������ ������ �������� � ���� ��������� ����������� ������
typedef struct Box_ {
	double x1, y1, z1, x2, y2, z2;
	int particles[100];
	short end;
} Box;

// ������ ������ ��� ����� ������
Box boxes_yz[16][16][64];


// ������ � �������� ������, ��� ������� ���� ��������� ������� �������
int particles_for_check[100];
int particles_for_check_count = 0;


void print_system_parameters() {
	long double etta = (PI * NP) / (6.0 * A * A * (L - 1.0));
	printf("\n\n| A    | %.16le |\n| L    | %.16le |\n| N    | %d                   |\n| etta | %.16le |\n", 2.0 * A, 2.0 * L, NP, etta);
}


void print_time_line() {
	for (int i = 0; i <= last; i++) {
		fprintf(stderr, "EVENT:");
		fprintf(stderr, "%d %d %d %.5le", i, time_queue[i].im, time_queue[i].jm, time_queue[i].t);
	}
}


double get_maximum_particle_time() {
	double t_max = -1.0e+20;

	for (int i = 0; i < NP; ++i) {
		if (particles[i].t > t_max) {
			t_max = particles[i].t;
		}
	}

	return t_max;
}


int check_particles() {
	double E = 0.0;
	double dt = 0.0;
	double t_global = get_maximum_particle_time();

	for (int i = 0; i < NP; i++) {
		particle p1 = particles[i];
		Box p1_box = boxes_yz[p1.y_box][p1.z_box][p1.x_box];

		dt = t_global - p1.t;
		p1.x += p1.vx * dt;
		p1.y += p1.vy * dt;
		p1.z += p1.vz * dt;

		E += p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;

		if ((p1.x_box > K2) || (p1.y_box > K) || (p1.z_box > K) ||
			(p1.x_box < 0) || (p1.y_box < 0) || (p1.z_box < 0)) {
			throw "Particle locates in incorrect cell.";
		}

		if (((abs(p1.x) - L) > 1.0e-14) ||
			((abs(p1.y) - A) > 1.0e-14) ||
			((abs(p1.z) - A) > 1.0e-14)) {
			printf(" \n Particle %d, %.15le, %.15le, %.15le, A = %.15le, %.15le \n ", i, p1.x, p1.y, p1.z, A, p1.z - A);
			throw "Particle is out of the system boundaries.";
		}

		if (((p1.x < p1_box.x1) && (abs(p1.x - p1_box.x1) > 1.0e-14)) ||
			((p1.x > p1_box.x2) && (abs(p1.x - p1_box.x2) > 1.0e-14)) ||
			((p1.y < p1_box.y1) && (abs(p1.y - p1_box.y1) > 1.0e-14)) ||
			((p1.y > p1_box.y2) && (abs(p1.y - p1_box.y2) > 1.0e-14)) ||
			((p1.z < p1_box.z1) && (abs(p1.z - p1_box.z1) > 1.0e-14)) ||
			((p1.z > p1_box.z2) && (abs(p1.z - p1_box.z2) > 1.0e-14))) {

			printf("Vilet za granicy %d \n", i);

			printf("Granizy:\n");
			printf("X : [%.15le ; %.15le]\n", p1_box.x1, p1_box.x2);
			printf("Y : [%.15le ; %.15le]\n", p1_box.y1, p1_box.y2);
			printf("Z : [%.15le ; %.15le]\n", p1_box.z1, p1_box.z2);
			printf("x, y, z: %.15le %.15le %.15le\n", p1.x, p1.y, p1.z);

			printf("p1.dt = %.15le, p.im = %d, p1.jm = %d \n", p1.t, time_queue[p1.ti].im, time_queue[p1.ti].jm);

			throw "Particle is out of the cell boundary.";
		}

		if ((p1.i_copy > -1) && (particles[i + NP].i_copy == -1)) {
			throw "Particle has incorrect image!";
		}

		bool w = false;
		for (int ty = 0; ty <= boxes_yz[p1.y_box][p1.z_box][p1.x_box].end; ++ty) {
			if (boxes_yz[p1.y_box][p1.z_box][p1.x_box].particles[ty] == i) w = true;
		}
		if (w == false) {
			for (int t = 0; t <= boxes_yz[p1.y_box][p1.z_box][p1.x_box].end; ++t) {
				printf(" %d ", boxes_yz[p1.y_box][p1.z_box][p1.x_box].particles[t]);
			}

			throw "Particle doesn't store in the cell.";
		}
		if (time_queue[p1.ti].im != i && time_queue[p1.ti].jm != i) {
			printf("\n >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ");
			printf("\n i = %d ; im = %d ; jm = %d ; ti = %d ", i, time_queue[p1.ti].im, time_queue[p1.ti].jm, p1.ti);
			throw "Particle has no correct link to the event.";
		}
	}

	if (abs(E - global_E) > 0.1e-8) {
		printf("\nENERGY was changed: \n E_seed = %.15le \n E_now= %.15le \n", global_E, E);
		throw "ENERGY was changed.";
	}

	for (int i = 0; i < last; i++) {
		if ((time_queue[i].im >= NP) && (particles[time_queue[i].im - NP].i_copy == -1)) {
			printf("\n Time tree event # %d ", i);
			printf("\n im, jm = %d %d \n", time_queue[i - 2].im, time_queue[i - 2].jm);
			printf("\n im, jm = %d %d \n", time_queue[i - 1].im, time_queue[i - 1].jm);
			printf("\n im, jm = %d %d \n", time_queue[i].im, time_queue[i].jm);
			printf("\n im, jm = %d %d \n", time_queue[i + 1].im, time_queue[i + 1].jm);
			printf("\n particle 34 event # %d \n", particles[34].ti);
			printf("\n last = %d \n", last);
			throw "Incorrect event!";
		}
		if ((time_queue[i].jm >= NP) && (particles[time_queue[i].jm - NP].i_copy == -1)) {
			throw "Incorrect event!";
		}
	}

	return 0;
}


/*
������� ������� �������� �� ������� ������� � ������ �������
���������:
i - �������, �� ������� ��������� ������� � ������ ������
t - ����� �� ����������� ������� �������
*/
int get_up(int i, double &t) {
	int j = i >> 1;
	while (i > 1) {
		if (i % 2 != 0 && t < time_queue[i - 1].t) {
			particles[time_queue[i - 1].im].ti = i;
			if (time_queue[i - 1].jm >= 0) particles[time_queue[i - 1].jm].ti = i;
			time_queue[i] = time_queue[i - 1];
			i--;
		}
		if (time_queue[j].t > t) {
			particles[time_queue[j].im].ti = i;
			if (time_queue[j].jm >= 0) particles[time_queue[j].jm].ti = i;
			time_queue[i] = time_queue[j];
			i = j;
			j >>= 1;
		}
		else return i;
	}
	return 1;
}


/*
������� ���������� ������� � ������� �������
���������:
e - ����� ������� (��. ��������� ��������� ���� Event)
*/
void add_event(int &i, int &j) {
	double t = particles[i].t + particles[i].dt;

	particles[i].ti = get_up(last, t);

	if (j >= 0) {
		particles[j].dt = particles[i].dt;
		particles[j].ti = particles[i].ti;
	}

	time_queue[particles[i].ti].im = i;
	time_queue[particles[i].ti].jm = j;
	time_queue[particles[i].ti].t = t;

	last++;
}


/*
������� �������� ������� �� ������� �������
���������:
i - ������� ���������� �������� � ������� �������
*/
void delete_event(int i) {
	int j = i << 1;
	while (j < last) {
		if (i % 2 == 0 && time_queue[i + 1].t < time_queue[j].t) {
			particles[time_queue[i + 1].im].ti = i;
			if (time_queue[i + 1].jm >= 0) particles[time_queue[i + 1].jm].ti = i;
			time_queue[i] = time_queue[i + 1];
			i++;
			j = i << 1;
		}
		else {
			particles[time_queue[j].im].ti = i;
			if (time_queue[j].jm >= 0) particles[time_queue[j].jm].ti = i;
			time_queue[i] = time_queue[j];
			i = j;
			j = i << 1;
		}
	}

	if (i < last - 1 && i % 2 == 0) {
		particles[time_queue[i + 1].im].ti = i;
		if (time_queue[i + 1].jm >= 0) particles[time_queue[i + 1].jm].ti = i;
		time_queue[i] = time_queue[i + 1];
		i++;
	}

	if (i < last - 1) {
		j = get_up(i, time_queue[last - 1].t);
		particles[time_queue[last - 1].im].ti = j;
		if (time_queue[last - 1].jm >= 0) particles[time_queue[last - 1].jm].ti = j;
		time_queue[j] = time_queue[last - 1];
	}

	last--;
}


void clear_particle_events(int &i) {
	// clear information about events for this virtual particle
	int f = particles[i].ti;

	if (f > 0) {
		int e = -100;
		int kim = time_queue[f].im;
		int kjm = time_queue[f].jm;

		if (time_queue[f].im == i) {
			if (time_queue[f].jm >= 0) {
				double dt = particles[kim].t - particles[kjm].t;  // ������� � ������� �������

				particles[kjm].x += particles[kjm].vx*dt;
				particles[kjm].y += particles[kjm].vy*dt;
				particles[kjm].z += particles[kjm].vz*dt;
				particles[kjm].t = particles[kim].t;
				particles[kjm].dt = (particles[kjm].dt - dt) / 1.1;

				delete_event(f);
				add_event(kjm, e);
			}
			else delete_event(f);
		}
		else if (time_queue[f].jm == i) {
			double dt = particles[kjm].t - particles[kim].t;  // ������� � ������� �������

			particles[kim].x += particles[kim].vx*dt;
			particles[kim].y += particles[kim].vy*dt;
			particles[kim].z += particles[kim].vz*dt;
			particles[kim].t = particles[kjm].t;
			particles[kim].dt = (particles[kim].dt - dt) / 1.1;

			delete_event(f);
			add_event(kim, e);
		}
		particles[i].ti = 0;
	}
}


/*
������� ������� ���������� ������� ��� �������
���������:
i - ����� �������
*/
void retime(int &i) {
	particle p1 = particles[i];
	Box p1_box = boxes_yz[p1.y_box][p1.z_box][p1.x_box];
	int jm;
	double dt, dt_min;

	clear_particle_events(i);

	if (p1.vx < 0.0) {
		dt_min = (p1_box.x1 - p1.x) / p1.vx;
		jm = -2;
		if (p1.x_box == 1) {
			dt_min = (p1_box.x1 + 1.0 - p1.x) / p1.vx;
			jm = -1;
		}
	}
	else {
		dt_min = (p1_box.x2 - p1.x) / p1.vx;
		jm = -4;
		if (p1.x_box == K2 - 1) {
			dt_min = (p1_box.x2 - 1.0 - p1.x) / p1.vx;
			jm = -1;
		}
	}

	if (p1.vy < 0.0) {
		dt = (p1_box.y1 - p1.y) / p1.vy;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -5;
		}
		if ((p1.y_box == 1) && (i < NP)) {
			dt = (p1_box.y1 + 1.0 - p1.y) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -15;
			}
		}
	}
	else {
		dt = (p1_box.y2 - p1.y) / p1.vy;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -6;
		}
		if ((p1.y_box == K - 1) && (i < NP)) {
			dt = (p1_box.y2 - 1.0 - p1.y) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -16;
			}
		}
	}

	if (p1.vz < 0.0) {
		dt = (p1_box.z1 - p1.z) / p1.vz;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -7;
		}
		if ((p1.z_box == 1) && (i < NP)) {
			dt = (p1_box.z1 + 1.0 - p1.z) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -17;
			}
		}
	}
	else {
		dt = (p1_box.z2 - p1.z) / p1.vz;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -8;
		}
		if ((p1.z_box == K - 1) && (i < NP)) {
			dt = (p1_box.z2 - 1.0 - p1.z) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -18;
			}
		}
	}

	double temp, dx, dy, dz, dvx, dvy, dvz, d, dv, bij;
	int s, n, r, q, w;

	for (r = p1.x_box - 1; r < p1.x_box + 2; ++r)
		for (q = p1.y_box - 1; q < p1.y_box + 2; ++q)
			for (w = p1.z_box - 1; w < p1.z_box + 2; ++w) {

				if (q == -1) {
					continue;
				}
				if (q == K + 1) {
					continue;
				}
				if (w == -1) {
					continue;
				}
				if (w == K + 1) {
					continue;
				}

				if (boxes_yz[q][w][r].end > 90) {
					printf("\n q = %d, w = %d, r = %d, end = %d \n", q, w, r, boxes_yz[q][w][r].end);

					for (s = 0; s < 32; ++s) {
						int u = boxes_yz[q][w][r].particles[s];
						printf("particle %d - %d %d %d \n", u, particles[u].x_box, particles[u].y_box, particles[u].z_box);
					}

					printf("\n s = %d , p1_box.end = %d\n", boxes_yz[q][w][r].end, p1_box.end);
					printf("\n %d %d %d \n", q, w, r);
					printf("stop");

					exit(0);

					for (s = 0; s < boxes_yz[q][w][r].end; ++s) {
						int u = boxes_yz[q][w][r].particles[s];
						printf("particle %d - %d %d %d \n", u, particles[u].x_box, particles[u].y_box, particles[u].z_box);
					}
					printf("stop");
				}

				for (s = 0; s <= boxes_yz[q][w][r].end; ++s) {
					n = boxes_yz[q][w][r].particles[s];

					if (n == p1.i_copy) continue;

					particle p = particles[n];

					temp = p1.t - p.t;

					dvx = p.vx - p1.vx;
					dvy = p.vy - p1.vy;
					dvz = p.vz - p1.vz;
					dx = p.x + p.vx * temp - p1.x;
					dy = p.y + p.vy * temp - p1.y;
					dz = p.z + p.vz * temp - p1.z;

					bij = dx * dvx + dy * dvy + dz*dvz;
					if (bij < 0.0) {
						dv = dvx * dvx + dvy * dvy + dvz*dvz;
						d = bij * bij + dv * (4.0 - dx * dx - dy * dy - dz * dz);

						if (d > 0.0) {
							dt = -(sqrt(d) + bij) / dv;

							/*
							��������, ��� ������� �������� ������������� �������:
							1. �����, ����������� � ������� ��� ����������� ���������� ����� ���������
							� ������ �������, ����� ����������� ���������� � �������� �� ������ ����
							����� ����� ���������.
							2. ������� n1 ������������ � �������� n2, �� ������� ����� ��� ������� n1,
							��� �������� �� ��������� ����� � �� ���������� ��� � ������ �������� n3,
							� ���������� ���� �������� ������� n1 ����� �������� � ����� ����������
							������ n1 � n2 ����� ���� ������������� (-1*10-14) �� �� �����������
							� ������� ��������� � 15�� �����.
							*/
							if (dt > -1.0e-12 && dt < 1.0e-15) dt = 0.0;

							temp += dt;

							if ((dt < dt_min) && (dt >= 0.0) &&
								((temp < p.dt) || ((abs(temp - p.dt) < 0.1E-15) &&
									(time_queue[p.ti].im == n) &&
									(time_queue[p.ti].jm == -100)))) {
								dt_min = dt;
								jm = n;
							}
						}
					}
				}
			}

	if (jm >= 0) {
		dt = p1.t - particles[jm].t;
		if (dt < 0.0 && dt > -1.0e-15)
			dt = 0.0;

		particles[jm].t = p1.t;
		particles[jm].dt = dt_min;
		particles[jm].x += particles[jm].vx * dt;
		particles[jm].y += particles[jm].vy * dt;
		particles[jm].z += particles[jm].vz * dt;

		clear_particle_events(jm);
	}

	particles[i].dt = dt_min;

	add_event(i, jm);

	if (dt_min < -1.0e-11) {
		printf("\n retime result: %d %d, %.16le\n ", i, jm, dt_min);
		printf("\n p1.x = %.16le, p1.y = %.16le, p1.z = %.16le \n", p1.x, p1.y, p1.z);
		printf("\n p1.x_box = %d, p1.y_box = %d, p1.z_box = %d", p1.x_box, p1.z_box, p1.y_box);
		printf("\n im = %d, jm = %d, dt = %.16le, A = %.16le", i, jm, dt_min, A);
		printf("\n p1.box.x = [%.16le; %.16le]", boxes_yz[p1.y_box][p1.z_box][p1.x_box].x1, boxes_yz[p1.y_box][p1.z_box][p1.x_box].x2);
		printf("\n p1.box.y = [%.16le; %.16le]", boxes_yz[p1.y_box][p1.z_box][p1.x_box].y1, boxes_yz[p1.y_box][p1.z_box][p1.x_box].y2);
		printf("\n p1.box.z = [%.16le; %.16le]", boxes_yz[p1.y_box][p1.z_box][p1.x_box].z1, boxes_yz[p1.y_box][p1.z_box][p1.x_box].z2);
		throw "dt < 0!";
	}
}


int search_collission_for_new_virtual_particle(int &i) {
	particle p1 = particles[i];
	particle p2;
	double dt, dx, dy, dz, dvx, bij, d, dv, dvy, dvz, fa, fb, fc, fD, sD, t1, t2;

	for (int j = 0; j < NP; j++) {

		if (j == i - NP) continue;

		p2 = particles[j];

		dt = p1.t - p2.t;

		// EN: if we have small time differences for these two particles
		// we shouldn't use it for collission because other virtual particle
		// can use the same particle in the same time.
		// RU: ���� � ������ ��� ������� �� ������� �� �� �� ������ ������������� ������������
		// ���� ������
		if (dt == 0.0) continue;

		p2.x += p2.vx*dt;
		p2.y += p2.vy*dt;
		p2.z += p2.vz*dt;

		fa = p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;
		fb = 2.0*p1.x*p1.vx + 2.0*p1.y*p1.vy + 2.0*p1.z*p1.vz - 2.0*p2.x*p1.vx - 2.0*p2.y*p1.vy - 2.0*p2.z*p1.vz;
		fc = p2.x*p2.x + p1.x*p1.x + p2.y*p2.y + p1.y*p1.y + p2.z*p2.z + p1.z*p1.z - 2.0*p1.x*p2.x - 2.0*p1.y*p2.y - 2.0*p1.z*p2.z - 4.0;

		fD = fb*fb - 4.0*fa*fc;
		if (fD > 0) {
			sD = sqrt(fD);

			t1 = (-fb - sD) / (2.0*fa);
			t2 = (-fb + sD) / (2.0*fa);

			dx = p2.x - p1.x - p1.vx*t1;
			dy = p2.y - p1.y - p1.vy*t1;
			dz = p2.z - p1.z - p1.vz*t1;

			dvx = p2.vx - p1.vx;
			dvy = p2.vy - p1.vy;
			dvz = p2.vz - p1.vz;

			bij = dx * dvx + dy * dvy + dz*dvz;
			if (bij < 0.0) {
				dv = dvx * dvx + dvy * dvy + dvz*dvz;
				d = bij * bij + dv * (4.0 - dx * dx - dy * dy - dz * dz);

				dt = -(sqrt(d) + bij) / dv;

				if (dt > 0.0) {
					particles[i].x += p1.vx*t1;
					particles[i].y += p1.vy*t1;
					particles[i].z += p1.vz*t1;

					particles[j].x = p2.x;
					particles[j].y = p2.y;
					particles[j].z = p2.z;

					dt = particles[i].t - particles[j].t;
					particles[j].t = particles[i].t;
					particles[j].dt -= dt;

					return j;
				}
			}
			else {
				dx = p2.x - p1.x - p1.vx*t2;
				dy = p2.y - p1.y - p1.vy*t2;
				dz = p2.z - p1.z - p1.vz*t2;

				bij = dx * dvx + dy * dvy + dz*dvz;
				if (bij < 0.0) {
					dv = dvx * dvx + dvy * dvy + dvz*dvz;
					d = bij * bij + dv * (4.0 - dx * dx - dy * dy - dz * dz);

					dt = -(sqrt(d) + bij) / dv;

					if (dt > 0.0) {
						particles[i].x += p1.vx*t2;
						particles[i].y += p1.vy*t2;
						particles[i].z += p1.vz*t2;

						particles[j].x = p2.x;
						particles[j].y = p2.y;
						particles[j].z = p2.z;

						dt = particles[i].t - particles[j].t;
						particles[j].t = particles[i].t;
						particles[j].dt -= dt;

						return j;
					}
				}
			}
		}
	}

	// EN: if we can't find collission we will return -1 and it will mean
	// than we should try to search collissions with other value of X for this virtual particle.
	// RU: ���� �� �� ����� ���������� ��� ����� ������, �� ���� �������� �
	// � ����������� ������ ���������� �����.
	return -1;

}

// ������ ������� ��������� ����� ��������� ����� ��� ������ ������
// �� ��� �. ���� ����� �� �������, ������� ���������� ����� �����
// � ����� �� ������, �������� ��� ���������
void find_place_for_particle(int &i) {

	double x, dy, dz, dx, d, dt, dvx, dvy, dvz, bij, dv, x_min, x_max = 0.0;
	double no_free_space_min[300];
	double no_free_space_max[300];

	int particle_on_the_line = -1;
	double particle_x_for_collission;

	bool include;
	int spaces = 2;
	no_free_space_min[1] = -L;
	no_free_space_max[1] = 1.0 - L + 1.0e-6;
	no_free_space_min[2] = L - 1.0 - 1.0e-6;
	no_free_space_max[2] = L;

	// BUG: ������������� �������� ���������� ������� �������� � 0 �� ����� �������
	// �.�. �� ��������� ���� ��� ����������� ������ ���������� �������, ��� �����������
	// ����������� ���� ��� ����� ������� ����� ������� ������ ����� ��������� ������.
	// ���������� ���������� ������ ���� ���������� ������ � �������� �������� ���� �� ���.
	for (int x_box = 0; x_box <= K2; ++x_box)
		for (int y_box = particles[i].y_box - 1; y_box < particles[i].y_box + 2; ++y_box) {
			for (int z_box = particles[i].z_box - 1; z_box < particles[i].z_box + 2; ++z_box) {

				if (y_box == -1) {
					continue;
				}
				if (y_box == K + 1) {
					continue;
				}
				if (z_box == -1) {
					continue;
				}
				if (z_box == K + 1) {
					continue;
				}

				for (int s = 0; s <= boxes_yz[y_box][z_box][x_box].end; ++s) {

					int n = boxes_yz[y_box][z_box][x_box].particles[s];

					/* we should ignore collisions of particles with itself: */
					if (n == i) continue;

					include = false;

					// calculate delta t between two particles.
					// If temp < 1.0e-14 we shouldn't use this particle for collission
					// because other virtual particle can use the same particle.
					double temp = particles[i].t - particles[n].t;

					dy = particles[n].y + particles[n].vy * temp - particles[i].y;
					dz = particles[n].z + particles[n].vz * temp - particles[i].z;
					d = 4.0 - dy*dy - dz*dz;
					if (d >= 0) {
						x_min = particles[n].x + particles[n].vx * temp - sqrt(d);
						x_max = particles[n].x + particles[n].vx * temp + sqrt(d);

						// search particle which can be used for collision
						if (((particle_on_the_line == -1) || (particle_on_the_line > NP)) &&
							(particles[i].vx * particles[n].vx < 0) && (temp > 1.0e-14)) {

							if (particles[i].vx < 0.0) x = x_max;
							else x = x_min;

							dx = particles[n].x + particles[n].vx * temp - x;
							dvx = particles[n].vx - particles[i].vx;
							dvy = particles[n].vy - particles[i].vy;
							dvz = particles[n].vz - particles[i].vz;

							bij = dx * dvx + dy * dvy + dz*dvz;
							if (bij < 0.0) {
								dv = dvx * dvx + dvy * dvy + dvz*dvz;
								d = bij * bij + dv * (4.0 - dx * dx - dy * dy - dz * dz);

								dt = -(sqrt(d) + bij) / dv;

								if ((dt > 0.0) && (abs(x) < L - 1.0)) {
									particle_on_the_line = n;
									particle_x_for_collission = x;
								}
							}
						}
						////////

						for (int m = 1; m <= spaces; ++m) {
							if (x_min <= no_free_space_min[m] &&
								x_max <= no_free_space_max[m] &&
								x_max >= no_free_space_min[m])
							{
								no_free_space_min[m] = x_min;
								include = true;
							}
							if (x_min <= no_free_space_min[m] &&
								x_max >= no_free_space_max[m])
							{
								no_free_space_min[m] = x_min;
								no_free_space_max[m] = x_max;
								include = true;
							}
							if (x_min <= no_free_space_max[m] &&
								x_min >= no_free_space_min[m] &&
								x_max >= no_free_space_max[m])
							{
								no_free_space_max[m] = x_max;
								include = true;
							}
						}
						if (include == false) {
							spaces += 1;
							no_free_space_min[spaces] = x_min;
							no_free_space_max[spaces] = x_max;
						}
					}
				}
			}
		}

	for (int r1 = 1; r1 < spaces; r1++)
		for (int r2 = r1 + 1; r2 <= spaces; r2++)
			if (no_free_space_min[r1] > no_free_space_min[r2])
			{
				double temp = no_free_space_min[r1];
				no_free_space_min[r1] = no_free_space_min[r2];
				no_free_space_min[r2] = temp;

				temp = no_free_space_max[r1];
				no_free_space_max[r1] = no_free_space_max[r2];
				no_free_space_max[r2] = temp;
			}

	/* we need to combine different free spaces if they overlap with each other */
	int r1 = 1;
	while (r1 < spaces)
	{
		// if two spaces overlap with each other or if we have no enough free space 
		// between two areas, we need to combine them
		if (no_free_space_max[r1] > no_free_space_min[r1 + 1])
		{
			if (no_free_space_max[r1] < no_free_space_max[r1 + 1])
				no_free_space_max[r1] = no_free_space_max[r1 + 1];
			for (int r2 = r1 + 1; r2 < spaces; r2++)
			{
				no_free_space_min[r2] = no_free_space_min[r2 + 1];
				no_free_space_max[r2] = no_free_space_max[r2 + 1];
			}
			spaces--;
		}
		else r1++;
	}

	if (spaces < 2) {    // if free space for virtual particles was not found:
		if (particle_on_the_line > -1) {

			particles[i].x = particle_x_for_collission;

			particle p = particles[particle_on_the_line];
			dt = particles[i].t - p.t; // time difference between two particles
									   // we need to move the second particle to
									   // the same time
			p.x += p.vx * dt;
			p.y += p.vy * dt;
			p.z += p.vz * dt;
			p.t = particles[i].t;
			p.dt -= dt;
			particles[particle_on_the_line] = p;

			return;
		}
		else {
			int r = -1, u = 0;
			
			r = search_collission_for_new_virtual_particle(i);

			// EN: randomize rand function
			// RU: ������������ ��������� ��������� �����
			srand(time(NULL));

			while (r == -1 && u < 1000)
			{
				// �������� ��������� ������� �� � ��� ������ � �������
				// ���, ����� ����� ��������� �� ������� ������ � ��������� �������
				double position_step = 2.0L * (L - 1.1) / double(RAND_MAX);
				particles[i].x = 1.1L + double(rand()) * position_step - L;

				// ���� ���������� � ������� ��������� ��� ����� ���������� � ��� ������
				// �.�. ����� ������������ ���������� �� ���������, � ������� ��������� ����� �������� ������
				// � ����� Y = const, Z = const, X in [1.1; L-1.1];
				r = search_collission_for_new_virtual_particle(i);

				u++;
			}

			// Stop the program if we can't find place for the virtual particle
			if (r == -1) {
				int parent_particle = i - NP;

				printf("\n A = %.15le, L = %.15le, u = %d \n ", A, L, u);

				printf("\n %d particle: x = %.15le, y = %.15le, z = %.15le \n", i, particles[i].x, particles[i].y, particles[i].z);
				printf("\n %d particle: vx = %.15le, vy = %.15le, vz = %.15le \n", i, particles[i].vx, particles[i].vy, particles[i].vz);
				printf("\n x_box = %d, y_box = %d, z_box = %d \n", particles[i].x_box, particles[i].y_box, particles[i].z_box);

				printf("\n %d particle: x = %.15le, y = %.15le, z = %.15le \n", parent_particle, particles[parent_particle].x, particles[parent_particle].y, particles[parent_particle].z);
				printf("\n %d particle: vx = %.15le, vy = %.15le, vz = %.15le \n", parent_particle, particles[parent_particle].vx, particles[parent_particle].vy, particles[parent_particle].vz);
				printf("\n x_box = %d, y_box = %d, z_box = %d \n", particles[parent_particle].x_box, particles[parent_particle].y_box, particles[parent_particle].z_box);

				throw "ERROR!!";
			}
		}
	}

	x_max = -L;
	x_min = L;
	double length = 2 * L, r;

	// select minimal free space and insert new particle in this space
	for (int m = 1; m < spaces; ++m) {
		// free space distance:
		r = abs(no_free_space_max[m] - no_free_space_min[m + 1]);
		if (r < length) {
			length = r;
			x_min = no_free_space_max[m];
			x_max = no_free_space_min[m + 1];
		}
	}

	particles[i].x = x_min + (x_max - x_min) / 2.0;
}


void destroy_virt_particle(int &i) {

	//FILE *history_file = fopen("history.txt", "a");
	//fprintf(history_file, "destroy virt_particle %d, x = %.15le, y = %.15le, z = %.15le\n", i, particles[i].x, particles[i].y, particles[i].z);
	//fclose(history_file);

	if (i >= NP) {
		printf("\n %d particle \n", i);
		printf(" i_copy = %d   particle icopy = %d \n", particles[i].i_copy, particles[i - NP].i_copy);
		printf(" x, y, z: %.5le %.5le %.5le", particles[i].x, particles[i].y, particles[i].z);

		throw "Error: can't destroy the real particle!";
	}
	if (particles[i].i_copy == -1) return;

	int new_i = i + NP;
	int x_box, y_box, z_box, box_i;

	x_box = particles[new_i].x_box;
	y_box = particles[new_i].y_box;
	z_box = particles[new_i].z_box;
	box_i = particles[new_i].box_i;
	Box p_box = boxes_yz[y_box][z_box][x_box];

	// clear information about this virtual particle in small cell
	if (p_box.particles[box_i] == new_i)
	{
		int j = p_box.particles[box_i] = p_box.particles[p_box.end];
		particles[j].box_i = box_i;
		--p_box.end;
		boxes_yz[y_box][z_box][x_box] = p_box;
	}

	particles[new_i].t = particles[i].t;
	clear_particle_events(new_i);

	particles[i].i_copy = -1;
	particles[new_i].i_copy = -1;
}


void change_with_virt_particles(int &im, int &jm) {
	int y_box, z_box;
	int f = im + NP;

	//FILE *history_file = fopen("history.txt", "a");
	//fprintf(history_file, "changing with virt_particle %d, x = %.15le, y = %.15le, z = %.15le\n", im, particles[im].x, particles[im].y, particles[im].z);
	//fclose(history_file);

	double dt = particles[im].t - particles[f].t;
	particles[im].x = particles[f].x + dt*particles[f].vx;
	particles[im].y = particles[f].y + dt*particles[f].vy;
	particles[im].z = particles[f].z + dt*particles[f].vz;

	// ���������� � � ����� ������
	y_box = particles[f].y_box;
	z_box = particles[f].z_box;

	if (y_box <= 0) y_box = 1;
	if (y_box >= K) y_box = K - 1;
	if (z_box <= 0) z_box = 1;
	if (z_box >= K) z_box = K - 1;

	particles[im].x_box = particles[f].x_box;
	particles[im].y_box = y_box;
	particles[im].z_box = z_box;

	// ��� ������ ������������ ������ �� ������� ������ ������� ����� �� ������� ������,
	// ����� �������� ���������� ������
	if (particles[im].y > A) particles[im].y = boxes_yz[y_box][z_box][particles[f].x_box].y2;
	else if (particles[im].y < -A) particles[im].y = boxes_yz[y_box][z_box][particles[f].x_box].y1;
	if (particles[im].z > A) particles[im].z = boxes_yz[y_box][z_box][particles[f].x_box].z2;
	else if (particles[im].z < -A) particles[im].z = boxes_yz[y_box][z_box][particles[f].x_box].z1;

	Box b = boxes_yz[y_box][z_box][particles[f].x_box];
	if (((particles[im].x < b.x1) && (abs(particles[im].x - b.x1) > 1.0e-14)) ||
		((particles[im].x > b.x2) && (abs(particles[im].x - b.x2) > 1.0e-14)) ||
		((particles[im].y < b.y1) && (abs(particles[im].y - b.y1) > 1.0e-14)) ||
		((particles[im].y > b.y2) && (abs(particles[im].y - b.y2) > 1.0e-14)) ||
		((particles[im].z < b.z1) && (abs(particles[im].z - b.z1) > 1.0e-14)) ||
		((particles[im].z > b.z2) && (abs(particles[im].z - b.z2) > 1.0e-14))) {
		printf("\n particle %d i_copy %d x = %.15le, y = %.15le, z = %.15le \n", im, particles[im].i_copy, particles[im].x, particles[im].y, particles[im].z);
		printf(" particle %d x = %.15le, y = %.15le, z = %.15le \n", f, particles[f].x, particles[f].y, particles[f].z);
		printf("box x [%.15le; %.15le] \n", b.x1, b.x2);
		printf("box y [%.15le; %.15le] \n", b.y1, b.y2);
		printf("box z [%.15le; %.15le] \n", b.z1, b.z2);

		printf("\n\n particle %d t = %.15le, particle %d t = %.15le \n", im, particles[im].t, f, particles[f].t);
		printf(" particle %d t+dt = %.15le, particle %d t+dt = %.15le \n", im, particles[im].t + particles[im].dt, f, particles[f].t + particles[f].dt);
		printf(" particle %d event: %d %d %.15le \n", f, time_queue[particles[f].ti].im, time_queue[particles[f].ti].jm, time_queue[particles[f].ti].t);

		throw "stop";
	}

	destroy_virt_particle(im);
}


void create_virt_particle(int &i, bool need_to_check=true) {
	double dt, dt_min, t_min, y, z, dy, dz;
	double kv = 1.0;
	double t01, t02;
	short x_box, y_box, z_box;
	int new_i = i + NP;

	destroy_virt_particle(i);

	//FILE *history_file = fopen("history.txt", "a");
	//fprintf(history_file, "trying to create virt_particle %d, x = %.15le, y = %.15le, z = %.15le\n", i, particles[i].x, particles[i].y, particles[i].z);
	//fclose(history_file);

	if ((need_to_check == false) ||
		(((particles[i].y <= 1.0 - A) && (particles[i].vy < 0.0)) ||
		 ((particles[i].y >= A - 1.0) && (particles[i].vy > 0.0)) ||
		 ((particles[i].z <= 1.0 - A) && (particles[i].vz < 0.0)) ||
		 ((particles[i].z >= A - 1.0) && (particles[i].vz > 0.0)))) {

		dt_min = 1.0e+20;

		//FILE *history_file = fopen("history.txt", "a");
		//fprintf(history_file, "creating virt_particle %d \n", i, particles[i].x, particles[i].y, particles[i].z);
		//fclose(history_file);

		y = A + particles[i].y;
		z = A + particles[i].z;
		dy = A + particles[i].y;
		dz = A + particles[i].z;

		if (particles[i].vy < 0)
			y = A - particles[i].y;
		if (particles[i].vz < 0)
			z = A - particles[i].z;

		if (particles[i].z > 0)
			dz = A - particles[i].z;
		if (particles[i].y > 0)
			dy = A - particles[i].y;

		if ((particles[i].y*particles[i].vy > 0) && (dy - 1.0 < 1.0e-14)) {
			dt = dy / abs(particles[i].vy);

			particles[new_i].vx = particles[i].vx;
			particles[new_i].vy = particles[i].vy;
			particles[new_i].vz = particles[i].vz;

			t01 = y / abs(particles[i].vy);
			t02 = z / abs(particles[i].vz);

			if (t01 < t02) {
				t_min = t01;

				particles[new_i].x = particles[i].x - particles[i].vx*t01 - particles[i].vx*dt;
				particles[new_i].y = particles[i].y - particles[i].vy*t01 - particles[i].vy*dt;
				particles[new_i].z = particles[i].z - particles[i].vz*t01 - particles[i].vz*dt;
			}
			else {
				t_min = t02;

				particles[new_i].x = particles[i].x - particles[i].vx*t02 - particles[new_i].vx*dt;
				particles[new_i].y = particles[i].y - particles[i].vy*t02 - particles[new_i].vy*dt;
				particles[new_i].z = particles[i].z - particles[i].vz*t02 - particles[new_i].vz*dt;
			}

			dt_min = dt;
		}
		if ((particles[i].z*particles[i].vz > 0) && (dz - 1.0 < 1.0e-14)) {
			dt = dz / abs(particles[i].vz);

			if (dt < dt_min) {

				particles[new_i].vx = particles[i].vx;
				particles[new_i].vy = particles[i].vy;
				particles[new_i].vz = particles[i].vz;

				t01 = y / abs(particles[i].vy);
				t02 = z / abs(particles[i].vz);

				if (t01 < t02) {
					t_min = t01;

					particles[new_i].x = particles[i].x - particles[i].vx*t01 - particles[i].vx*dt;
					particles[new_i].y = particles[i].y - particles[i].vy*t01 - particles[i].vy*dt;
					particles[new_i].z = particles[i].z - particles[i].vz*t01 - particles[i].vz*dt;
				}
				else {
					t_min = t02;

					particles[new_i].x = particles[i].x - particles[i].vx*t02 - particles[new_i].vx*dt;
					particles[new_i].y = particles[i].y - particles[i].vy*t02 - particles[new_i].vy*dt;
					particles[new_i].z = particles[i].z - particles[i].vz*t02 - particles[new_i].vz*dt;
				}

				dt_min = dt;
			}
		}

		// Find correct X coordinate for new particle
		if (abs(particles[new_i].x) > L - 1.0) {
			int f = 1;
			double LL = L - 1.0;
			double x = particles[new_i].x;

			if (particles[new_i].x < -LL) f = -1;
			x = x - f*LL;
			int m = int(x / (2.0 * LL));
			x = x - 2.0 * m * LL;
			int x_wall = 1;
			if (m % 2 != 0)
				x_wall = -1;
			else
				particles[new_i].vx = -particles[new_i].vx;
			particles[new_i].x = x_wall * (f * LL - x);
		}

		x_box = short((L + dL + particles[new_i].x) / dL);
		y_box = short((A + dA + particles[new_i].y) / dA);
		z_box = short((A + dA + particles[new_i].z) / dA);

		if (y_box < 0) y_box = 0;
		if (z_box < 0) z_box = 0;
		if (y_box > K) y_box = K;
		if (z_box > K) z_box = K;

		// ������������ ������� ������, � ������� �������� �������, �.�. ��� ������� double
		// �� ����� ����������� � 1.0e-13, �� �� ���� ������ ����� ���� ���������� �������
		if ((boxes_yz[y_box][z_box][x_box].x1 > particles[new_i].x) && (x_box > 0)) x_box--;
		if ((boxes_yz[y_box][z_box][x_box].x2 < particles[new_i].x) && (x_box < K2)) x_box++;
		if ((boxes_yz[y_box][z_box][x_box].y1 > particles[new_i].y) && (y_box > 0)) y_box--;
		if ((boxes_yz[y_box][z_box][x_box].y2 < particles[new_i].y) && (y_box < K)) y_box++;
		if ((boxes_yz[y_box][z_box][x_box].z1 > particles[new_i].z) && (z_box > 0)) z_box--;
		if ((boxes_yz[y_box][z_box][x_box].z2 < particles[new_i].z) && (z_box < K)) z_box++;

		particles[new_i].x_box = x_box;
		particles[new_i].y_box = y_box;
		particles[new_i].z_box = z_box;

		particles[new_i].t = particles[i].t;

		// �������� ���� ��� ����� ����� �� ������������ � ��� ������������� ���������
		bool search = false;
		for (short r = x_box - 1; r < x_box + 2; ++r) {
			for (short q = y_box - 1; q < y_box + 2; ++q) {
				for (short w = z_box - 1; w < z_box + 2; ++w) {

					if (q == -1) {
						continue;
					}
					if (q == K + 1) {
						continue;
					}
					if (w == -1) {
						continue;
					}
					if (w == K + 1) {
						continue;
					}

					for (int s = 0; s <= boxes_yz[q][w][r].end; ++s) {
						int n = boxes_yz[q][w][r].particles[s];
						double temp = particles[new_i].t - particles[n].t;
						double dx = particles[n].x + particles[n].vx * temp - particles[new_i].x;
						double dy = particles[n].y + particles[n].vy * temp - particles[new_i].y;
						double dz = particles[n].z + particles[n].vz * temp - particles[new_i].z;

						if (4.0 > dx*dx + dy*dy + dz*dz) {
							search = true;
						}
					}
				}
			}
		}

		if (search == true) {

			find_place_for_particle(new_i);

			/* update information about cell for new virtual particle */
			x_box = short((L + dL + particles[new_i].x) / dL);
			y_box = short((A + dA + particles[new_i].y) / dA);
			z_box = short((A + dA + particles[new_i].z) / dA);

			if (y_box < 0) y_box = 0;
			if (z_box < 0) z_box = 0;
			if (y_box > K) y_box = K;
			if (z_box > K) z_box = K;

			// ������������ ������� ������, � ������� �������� �������, �.�. ��� ������� double
			// �� ����� ����������� � 1.0e-13, �� �� ���� ������ ����� ���� ���������� �������
			if ((boxes_yz[y_box][z_box][x_box].x1 > particles[new_i].x) && (x_box > 0)) x_box--;
			if ((boxes_yz[y_box][z_box][x_box].x2 < particles[new_i].x) && (x_box < K2)) x_box++;
			if ((boxes_yz[y_box][z_box][x_box].y1 > particles[new_i].y) && (y_box > 0)) y_box--;
			if ((boxes_yz[y_box][z_box][x_box].y2 < particles[new_i].y) && (y_box < K)) y_box++;
			if ((boxes_yz[y_box][z_box][x_box].z1 > particles[new_i].z) && (z_box > 0)) z_box--;
			if ((boxes_yz[y_box][z_box][x_box].z2 < particles[new_i].z) && (z_box < K)) z_box++;

			particles[new_i].x_box = x_box;
			particles[new_i].y_box = y_box;
			particles[new_i].z_box = z_box;
		}

		short end = particles[new_i].box_i = ++boxes_yz[y_box][z_box][x_box].end;
		boxes_yz[y_box][z_box][x_box].particles[end] = new_i;

		particles[new_i].i_copy = i;
		particles[i].i_copy = new_i;
	}
}


void load_information_about_one_particle(FILE *loading_file) {
	double a1, a2, a3;
	int i, i_copy, x_box, y_box, z_box, end;

	fscanf(loading_file, "%d %d\n", &i, &i_copy);
	particles[i].i_copy = i_copy;
	fscanf(loading_file, "%le %le %le\n", &a1, &a2, &a3);
	particles[i].x = a1 - L;
	particles[i].y = a2 - A;
	particles[i].z = a3 - A;

	x_box = short((a1 + dL) / dL);
	y_box = short((a2 + dA) / dA);
	z_box = short((a3 + dA) / dA);

	if (y_box < 0) y_box = 0;
	if (z_box < 0) z_box = 0;
	if (y_box > K) y_box = K;
	if (z_box > K) z_box = K;

	// ������������ ������� ������, � ������� �������� �������, �.�. ��� ������� double
	// �� ����� ����������� � 1.0e-13, �� �� ���� ������ ����� ���� ���������� �������
	if ((boxes_yz[y_box][z_box][x_box].x1 > particles[i].x) && (x_box > 0)) x_box--;
	if ((boxes_yz[y_box][z_box][x_box].x2 < particles[i].x) && (x_box < K2)) x_box++;
	if ((boxes_yz[y_box][z_box][x_box].y1 > particles[i].y) && (y_box > 0)) y_box--;
	if ((boxes_yz[y_box][z_box][x_box].y2 < particles[i].y) && (y_box < K)) y_box++;
	if ((boxes_yz[y_box][z_box][x_box].z1 > particles[i].z) && (z_box > 0)) z_box--;
	if ((boxes_yz[y_box][z_box][x_box].z2 < particles[i].z) && (z_box < K)) z_box++;

	particles[i].x_box = x_box;
	particles[i].y_box = y_box;
	particles[i].z_box = z_box;

	if (i < NP) {
		if (boxes_yz[y_box][z_box][x_box].x1 <= particles[i].x &&
			boxes_yz[y_box][z_box][x_box].x2 >= particles[i].x)
			particles[i].x_box = x_box;
		else {
			printf("Particle locates in incorrect place %d", i);
			printf("\n a1 = %.15le, x = %.15le, L = %.15le, dL = %.15le \n", a1, particles[i].x, L, dL);
			printf("\n x_box = %d, BOX x: [%.15le; %.15le]", x_box, boxes_yz[y_box][z_box][x_box].x1, boxes_yz[y_box][z_box][x_box].x2);
			printf("\n particle x: %.15le", particles[i].x);
			throw "Particle locates in incorrect place";
		}

		if (boxes_yz[y_box][z_box][x_box].y1 <= particles[i].y &&
			boxes_yz[y_box][z_box][x_box].y2 >= particles[i].y)
			particles[i].y_box = y_box;
		else {
			printf("Particle locates in incorrect place %d", i);
			printf("\n BOX y: [%.15le; %.15le]", boxes_yz[y_box][z_box][x_box].y1, boxes_yz[y_box][z_box][x_box].y2);
			printf("\n particle y: %.15le", particles[i].y);
			throw "Particle locates in incorrect place";
		}

		if (boxes_yz[y_box][z_box][x_box].z1 <= particles[i].z &&
			boxes_yz[y_box][z_box][x_box].z2 >= particles[i].z)
			particles[i].z_box = z_box;
		else {
			printf("\n Particle %d locates in incorrect place", i);
			printf("\n BOX z: [%.15le; %.15le]", boxes_yz[y_box][z_box][x_box].z1, boxes_yz[y_box][z_box][x_box].z2);
			printf("\n particle Z: %.15le", particles[i].z);
			throw "Particle locates in incorrect place";
		}
	}

	fscanf(loading_file, "%le %le %le\n", &a1, &a2, &a3);

	particles[i].vx = a1;
	particles[i].vy = a2;
	particles[i].vz = a3;
	particles[i].t = 0.0;
	particles[i].dt = 0.0;
	particles[i].ti = 1;
	end = particles[i].box_i = ++boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].end;
	boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].particles[particles[i].box_i] = i;

	global_E += a1*a1 + a2*a2 + a3*a3;


	if (end > 7) {
		printf("Error: Too many particles in one box");
		printf("\n x_box, y_box, z_box = %d %d %d , end = %d \n", x_box, y_box, z_box, end);
		throw "Error: Too many particles in one box";
	}

	for (short t = 0; t < particles[i].box_i; ++t)
		for (short d = t + 1; d <= particles[i].box_i; ++d)
			if (boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].particles[t] == boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].particles[d])
			{
				printf("\n Duplicated particles in one box: %d \n", i);

				printf("a");

				for (int j = 0; j < boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].end; ++j)
					printf("%d ", boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].particles[j]);

				throw "Duplicated particles in one box";
			}
}


/*
������� �������� ������ � ������� �� �����.
����������� ��������� ������, ���������� � �������� ���� ������
*/
void load_seed(std::string file_name) {
	double a1, x, y, z;
	int i, count_of_virt_particles;

	FILE *loading_file;
	try {
		loading_file = fopen(file_name.c_str(), "r");
	}
	catch (...) {
		printf("\n ERROR: Can't read file '%s' - is it exist? \n", file_name.c_str());
		exit(1);
	}

	fscanf(loading_file, "%i\n", &i);
	NP = i;
	fscanf(loading_file, "%i\n", &i);
	count_of_virt_particles = i;
	fscanf(loading_file, "%le\n", &a1);
	A = a1 / 2.0;
	A2 = a1;
	fscanf(loading_file, "%le\n", &a1);
	L = a1 / 2.0;

	// EN: search for the appropriate values for dA, dL, K, K2
	// RU: ���� ���������� �������� ��� dA, dL, K, K2
	dA = 2.5;
	dL = 2.5;
	K = short(A2 / dA) + 1;
	K2 = short(2.0 * L / dL) + 1;
	dA = A2 / (K - 1);
	dL = 2.0*L / (K2 - 1);

	y = z = -A - dA;
	x = -L - dL;
	for (int i = 0; i <= K; i++) {
		for (int j = 0; j <= K; j++) {
			for (int w = 0; w <= K2; w++) {
				boxes_yz[i][j][w].x1 = x;
				boxes_yz[i][j][w].x2 = x + dL;
				boxes_yz[i][j][w].y1 = y;
				boxes_yz[i][j][w].y2 = y + dA;
				boxes_yz[i][j][w].z1 = z;
				boxes_yz[i][j][w].z2 = z + dA;
				boxes_yz[i][j][w].end = -1;

				x += dL;
			}
			x = -L - dL;
			z += dA;
		}
		y += dA;
		z = -A - dA;
	}

	for (int i = 0; i < NP * 2; i++) {
		particles[i].x_box = 0;
		particles[i].y_box = 0;
		particles[i].z_box = 0;
		particles[i].box_i = 0;
		particles[i].i_copy = -1;
		particles[i].ti = 0;
		particles[i].dt = 0.0;
	}

	global_E = 0.0;
	for (int i = 0; i < NP + count_of_virt_particles; i++) {
		load_information_about_one_particle(loading_file);
	}

	fclose(loading_file);

	last = 1;
	time_queue[0].t = 0.0;
	time_queue[0].im = -1;
	for (int i = 1; i < 16384; ++i) time_queue[i].t = 1.0E+20;

	for (int i = 0; i < NP; ++i) {
		retime(i);
		if (particles[i].i_copy > 0) retime(particles[i].i_copy);
	}

	printf("\n System was successfully loaded from file '%s' \n", file_name.c_str());
}


/*
������� ���������� ��������� ������� � ��������� ����.
����� ����������� ������������ ������������� ������ �� �������.
���������:
file_name - ��� �����, � ������� ����� �������� ����������
*/
void save(std::string file_name) {
	int images[7000];
	int count_of_images = 0;
	double x, y, z, dt;
	double t_global = get_maximum_particle_time();

	// RU: �� ������ ��������� ���������� �� ������� ������ �����
	// �� ������������� ������ ��� ������ �������� ������� �� ������������ ���������.
	// EN: we need to save all virtual particles (images) because we
	// don't want to create all virtual particles during the load_seed().
	for (int i = 0; i < NP; ++i) {
		if (particles[i].i_copy >= NP) {
			images[count_of_images] = particles[i].i_copy;
			count_of_images += 1;
		}
	}

	FILE *save_file = fopen(file_name.c_str(), "w+");

	// RU: ��������� ���������� � ���������� ������ � �������� �������
	// EN: save information about count of particles and size of the system
	fprintf(save_file, "%d\n", NP);
	fprintf(save_file, "%d\n", count_of_images);
	fprintf(save_file, "%.15le\n", A * 2.0);
	fprintf(save_file, "%.15le\n", L * 2.0);

	// RU: ��������� ���������� ���� ������ � �� ��������: x, y, z, vx, vy, vz
	// EN: we need to save all coordinates of particles: x, y, z, vx, vy, vz
	for (int i = 0; i < NP; ++i) {
		fprintf(save_file, "%d %d\n", i, particles[i].i_copy);

		dt = t_global - particles[i].t;

		x = L + particles[i].x + particles[i].vx * dt;
		y = A + particles[i].y + particles[i].vy * dt;
		z = A + particles[i].z + particles[i].vz * dt;
		fprintf(save_file, "%.15le %.15le %.15le\n", x, y, z);
		fprintf(save_file, "%.15le %.15le %.15le\n", particles[i].vx, particles[i].vy, particles[i].vz);
	}

	// RU: �� ��� �� ��������� ���������� ���� ����������� ������ (�������)
	// EN: we also need to save all coordinates of virtual particles (images)
	for (int i = 0; i < count_of_images; ++i) {
		int m = images[i];
		fprintf(save_file, "%d %d\n", m, particles[m].i_copy);

		dt = t_global - particles[m].t;
		x = L + particles[m].x + particles[m].vx * dt;
		y = A + particles[m].y + particles[m].vy * dt;
		z = A + particles[m].z + particles[m].vz * dt;
		fprintf(save_file, "%.15le %.15le %.15le\n", x, y, z);
		fprintf(save_file, "%.15le %.15le %.15le\n", particles[m].vx, particles[m].vy, particles[m].vz);
	}

	fclose(save_file);
}


// ��������� ������ ������ �������.
// �������� � �������� ����������
// NN - ����� ������ � �����
// etta - ��������� ���������
void new_seed(int NN, double etta) {
	int KZ = 2;

	double axy, axz, v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z;

	double rb = 2.0 * sqrt(2.0); // ���������� ����� ����� ������
	double rbp = sqrt(2.0);
	double A0 = rb * (NN - 1.0) + 2.0; // ������� ���������� ���������� ������
	double L0 = rb * (2.0 * NN - 1.0) + 2.0; //

	double mL0 = L0;

	double Betta = 0.0;

	double XC, YC, ZC; // ���������� ������ ������
	int NA; // ���������� ��������� ������
	double ax, ay, az, vx, vy, vz;

	NP = 2.0 * NN * (NN * NN + (NN - 1.0) * (NN - 1.0)); //  ��������� ����� ������.
	NP = NP + (2.0 * NN - 1.0) * 2.0 * (NN - 1.0) * NN; //

	double etta0 = (4.0 * PI * NP) / (3.0 * A0 * A0 * (L0 - 2.0)); // ������� �������������� ���������.

	double Alpha = etta0 / etta;

	double sk = L0 / (2.0 * Alpha * (L0 - 2.0));

	Betta = exp((1.0 / 3.0) * log(L0 / (2.0 * Alpha * (L0 - 2.0)) + sk));
	//Betta = Betta + exp( (1.0 / 3.0) * log( -L0 / (2.0 * Alpha * (L0-2.0)) - sk) );

	Betta = 1.0 / Betta;

	XC = L0 / 2.0; //
	YC = A0 / 2.0; //  ���������� ������ ������.
	ZC = A0 / 2.0; //

	L0 = L0*KZ;
	L = L0 * Betta;  // ��������� Betta - ��� ����������� ����������
	A = A0 * Betta;  // �� ���� ���������� ��� ���������� � ��������� ������

	NA = 0;

	srand(time(NULL));

	for (int i = 0; i < 2 * NN; i++)
		for (int j = 0; j < NN; j++)
			for (int k = 0; k < NN / 2; k++) {
				vx = double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)) - double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)); //  ������ ��������.
				vy = double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)) - double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)); //  ��� ����� ��������� ��� ���������� ������ �����(� ��������� �����������)
				vz = double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)) - double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)); //  ����� ��� ����� ���� ������� � ������ �������� �������.

				for (int ii = -1; ii < 2; ii += 2) {
					if ((i == 0) && (ii > 0)) continue;

					ax = XC + i * ii*rbp + 1.0e-7;

					for (int jj = -1; jj < 2; jj += 2) {
						if ((j == 0) && (jj > 0)) continue;

						ay = YC + j * jj*rbp + 1.0e-7;

						for (int kk = -1; kk < 2; kk += 2) {
							if (((i % 2 == 0) && (j % 2 == 0)) || ((i % 2 != 0) && (j % 2 != 0)))
								az = ZC + kk * (rbp + k * rb) + 1.0e-7;
							else {
								if ((k == 0) && (kk > 0)) continue;
								az = ZC + kk * k*rb + 1.0e-7;
							}

							particles[NA].x = ax; // ��������� ���������� ��� ����� ��������
							particles[NA].y = ay;
							particles[NA].z = az;

							if ((j == 0) && (fabs(ZC - az) < 0.1E-15L) && ((i / 2) * 2 != i))
								if (((i + 1) / 4) * 4 != (i + 1)) {
									axy = vy;
									axz = vz;
									if (ii > 0) {
										vy = -vy;
										vz = -vz;
									}
									vy = ii*vy;
									vz = ii*vz;
								}
								else {
									vy = axy;
									vz = axz;
									vy = ii*vy;
									vz = ii*vz;
								}
								if (i == 0) {
									if ((j == 0) && (kk < 0)) {
										switch (k) {
										case 0:
											v0x = vx;
											v0y = vy;
											v0z = vz;
											break;
										case 1:
											v1x = vx;
											v1y = vy;
											v1z = vz;
											break;
										case 2:
											v2x = vx;
											v2y = vy;
											v2z = vz;
											break;
										case 3:
											v3x = vx;
											v3y = vy;
											v3z = vz;
											break;
										}
									}
									else {
										switch (k) {
										case 0:
											vx = -v0x;
											vy = -v0y;
											vz = -v0z;
											break;
										case 1:
											vx = -v1x;
											vy = -v1y;
											vz = -v1z;
											break;
										case 2:
											vx = -v2x;
											vy = -v2y;
											vz = -v2z;
											break;
										case 3:
											vx = -v3x;
											vy = -v3y;
											vz = -v3z;
											break;
										}
									}

									if ((k == 0) && ((j / 2) * 2 != j)) {
										switch (j) {
										case 1:
											vx = jj*v0x;
											vy = jj*v0y;
											vz = jj*v0z;
											break;
										case 3:
											vx = jj*v1x;
											vy = jj*v1y;
											vz = jj*v1z;
											break;
										case 5:
											vx = jj*v2x;
											vy = jj*v2y;
											vz = jj*v2z;
											break;
										case 7:
											vx = jj*v3x;
											vy = jj*v3y;
											vz = jj*v3z;
											break;
										}
									}
								}

								particles[NA].vx = vx * ii * jj * kk;
								particles[NA].vy = vy * ii * jj * kk;
								particles[NA].vz = vz * ii * jj * kk;

								NA++;
						}
					}
				}
			}

	// �������� ������� ������� ���, ������� �����������.
	for (int jj = 1; jj < KZ; jj++) {
		for (int ii = 0; ii < NP; ii++) {
			particles[ii + jj * NP].x = particles[ii].x + jj*mL0;
			particles[ii + jj * NP].y = particles[ii].y;
			particles[ii + jj * NP].z = particles[ii].z;

			particles[ii + jj * NP].vx = particles[ii].vx;
			particles[ii + jj * NP].vy = particles[ii].vy;
			particles[ii + jj * NP].vz = particles[ii].vz;
		}
	}

	NP = NP*KZ;

	// "���������" ������� �� ����������� ���������.
	for (int ii = 0; ii < NP; ii++) {
		particles[ii].x = particles[ii].x*Betta;   // ��������� Betta - ��� ����������� ����������
		particles[ii].y = particles[ii].y*Betta;   // �� ���� ���������� ��� ���������� � ��������� ������
		particles[ii].z = particles[ii].z*Betta;
	}

	double Lx = 0.0;
	for (int i = 0; i < NP; i++) {
		Lx += (particles[i].y + A)*particles[i].vz - (particles[i].z + A)*particles[i].vy;
	}

	printf("\n Lx = %.15le\n", Lx);

	A2 = A;
	dA = 2.5;
	dL = 2.5;
	K = short(A / dA) + 1;
	K2 = short(L / dL) + 1;
	dA = A / (K - 1);
	dL = L / (K2 - 1);

	A = A / 2.0;
	L = L / 2.0;
	for (int i = 0; i < NP; i++) {
		particles[i].t = 0.0;
		particles[i].dt = 0.0;
		particles[i].ti = 0;

		particles[i].x -= L;
		particles[i].y -= A;
		particles[i].z -= A;

		particles[i].i_copy = -1;
	}

	L = ((PI * NP) / etta) / (6.0 * A * A) + 1.0;  // �������� ��� ������� ������� �������� ��������� ���������

												 // ��������� � ��������� ������� �� ����� ����� ����������������
												 // ��� ����������� ����������
	save("new.txt");
	load_seed("new.txt");
}


/*
������� ��������� ��������� ������ � ������������ � �����������
� ������� ��������.
���������:
im - ����� �������
jm - ����� ������ ������� ��� ����� �������
*/
bool reform(int &im, int &jm) {
	particle p1 = particles[im];
	double q1, q2, z, dx, dy, dz;
	bool need_create_virt_particle = false;

	if (jm >= 0) {
		particle p2 = particles[jm];

		double dt = p2.dt;

		p2.x += p2.vx * p2.dt;
		p2.y += p2.vy * p2.dt;
		p2.z += p2.vz * p2.dt;
		p2.t = p1.t;
		p2.dt = 0.0;
		dx = p1.x - p2.x;
		dy = p1.y - p2.y;
		dz = p1.z - p2.z;

		if (im >= NP) {
			int m = im - NP;

			p1.vy = particles[m].vy;
			p1.vz = particles[m].vz;

			if (p1.vx * particles[m].vx < 0) {
				p1.vx = -particles[m].vx;
			}
			else {
				p1.vx = particles[m].vx;
			}
		}
		if (jm >= NP) {
			int m = jm - NP;

			p2.vy = particles[m].vy;
			p2.vz = particles[m].vz;

			if (p2.vx * particles[m].vx < 0) {
				p2.vx = -particles[m].vx;
			}
			else {
				p2.vx = particles[m].vx;
			}
		}

		q1 = (dx * p1.vx + dy * p1.vy + dz * p1.vz) / 4.0;
		q2 = (dx * p2.vx + dy * p2.vy + dz * p2.vz) / 4.0;
		z = q2 - q1;
		p1.vx += dx*z;
		p1.vy += dy*z;
		p1.vz += dz*z;
		z = q1 - q2;
		p2.vx += dx*z;
		p2.vy += dy*z;
		p2.vz += dz*z;

		int e1 = im;
		if (im >= NP) {
			e1 = im - NP;
			particles[e1].vx = p1.vx;
			particles[e1].vy = p1.vy;
			particles[e1].vz = p1.vz;
		}
		else
			particles[im] = p1;

		int e2 = jm;
		if (jm >= NP) {
			e2 = jm - NP;

			double delta = p2.t - particles[e2].t;
			particles[e2].x += particles[e2].vx * delta;
			particles[e2].y += particles[e2].vy * delta;
			particles[e2].z += particles[e2].vz * delta;
			particles[e2].t = p2.t;

			particles[e2].vx = p2.vx;
			particles[e2].vy = p2.vy;
			particles[e2].vz = p2.vz;
		}
		else
			particles[jm] = p2;

		if (particles[im].i_copy > 0) {
			int m = im;
			if (particles[im].i_copy < NP)
				m = m - NP;
			destroy_virt_particle(m);
			p1.i_copy = -1;
		}
		if (particles[jm].i_copy > 0) {
			int m = jm;
			if (particles[jm].i_copy < NP)
				m = m - NP;
			destroy_virt_particle(m);
			p2.i_copy = -1;
		}

		create_virt_particle(e1);
		create_virt_particle(e2);
	}
	else
		if (jm == -1) {
			p1.vx = -p1.vx;
			particles[im] = p1;
		}
		else {
			if (jm != -100) {
				short end = boxes_yz[p1.y_box][p1.z_box][p1.x_box].end;
				Box p1_box = boxes_yz[p1.y_box][p1.z_box][p1.x_box];

				// ������� ������� �� ������ � ������� ��� ����������,
				// ������ �� ��������� ������� �� ����� ������ ������ � ������ �
				// ��������� ����� ������ � ������ �� 1.
				p1_box.particles[p1.box_i] = p1_box.particles[end];
				particles[p1_box.particles[end]].box_i = p1.box_i;
				--p1_box.end;
				// save
				boxes_yz[p1.y_box][p1.z_box][p1.x_box] = p1_box;

				if (jm == -2) {
					--p1.x_box;
					p1.x = boxes_yz[p1.y_box][p1.z_box][p1.x_box].x2;
				}
				if (jm == -4) {
					++p1.x_box;
					p1.x = boxes_yz[p1.y_box][p1.z_box][p1.x_box].x1;
				}
				if (jm == -5) {
					--p1.y_box;
					p1.y = boxes_yz[p1.y_box][p1.z_box][p1.x_box].y2;
					if ((p1.y_box == 0) && (im < NP)) {
						change_with_virt_particles(im, jm);
						p1 = particles[im];
						need_create_virt_particle = true;
					}
				}
				if (jm == -6) {
					++p1.y_box;
					p1.y = boxes_yz[p1.y_box][p1.z_box][p1.x_box].y1;
					if ((p1.y_box == K) && (im < NP)) {
						change_with_virt_particles(im, jm);
						p1 = particles[im];
						need_create_virt_particle = true;
					}
				}
				if (jm == -7) {
					--p1.z_box;
					p1.z = boxes_yz[p1.y_box][p1.z_box][p1.x_box].z2;
					if ((p1.z_box == 0) && (im < NP)) {
						change_with_virt_particles(im, jm);
						p1 = particles[im];
						need_create_virt_particle = true;
					}
				}
				if (jm == -8) {
					++p1.z_box;
					p1.z = boxes_yz[p1.y_box][p1.z_box][p1.x_box].z1;
					if ((p1.z_box == K) && (im < NP)) {
						change_with_virt_particles(im, jm);
						p1 = particles[im];
						need_create_virt_particle = true;
					}
				}
				if (jm < -10) {

					//if (jm == -15) particles[im].y = boxes_yz[p1.y_box][p1.z_box][p1.x_box].y1 + 1.0;
					//if (jm == -16) particles[im].y = boxes_yz[p1.y_box][p1.z_box][p1.x_box].y2 - 1.0;
					//if (jm == -17) particles[im].z = boxes_yz[p1.y_box][p1.z_box][p1.x_box].z1 + 1.0;
					//if (jm == -18) particles[im].z = boxes_yz[p1.y_box][p1.z_box][p1.x_box].z2 - 1.0;

					// ������ ����� �� �������� ��������� �������:
					create_virt_particle(im, false);
					p1 = particles[im];
				}

				p1.box_i = ++boxes_yz[p1.y_box][p1.z_box][p1.x_box].end;
				boxes_yz[p1.y_box][p1.z_box][p1.x_box].particles[p1.box_i] = im;
				particles[im] = p1;

				//printf("\n 4 end = %d \n", boxes_yz[p1.y_box][p1.z_box][p1.x_box].end);
			}
		}

		if ((p1.i_copy > -1) && (jm == -1))
		{
			int e = im;
			if (im >= NP) {
				e = im - NP;
				particles[e].vx = -particles[e].vx;
			}
			clear_particle_events(e);
			create_virt_particle(e);
		}

		if (need_create_virt_particle == true) {
			int e = im;
			if (im >= NP)
				e = im - NP;
			clear_particle_events(e);
			create_virt_particle(e);
		}

		return need_create_virt_particle;
}


/*
������� "���", �������� ���� ���������
*/
void step() {
	particle p1;
	int i, im, jm;
	double time = 0.0;
	bool need_virt_particle_retime;

	COLL_COUNT = 0;
	jm = 0;

	while (COLL_COUNT < NP / 2 || jm != -100) {
		im = time_queue[1].im;
		jm = time_queue[1].jm;

		//check_particles();
		//FILE *history_file = fopen("history.txt", "a");
		//fprintf(history_file, "im = %d, jm = %d, p[im].icopy = %d, p[jm].i_copy = %d, dt = %.16le\n", im, particles[im].i_copy, jm, particles[jm].i_copy, particles[im].dt);
		//fclose(history_file);

		delete_event(1);

		p1 = particles[im];

		p1.ti = -1;
		if (jm >= 0) particles[jm].ti = -1;

		if (jm == -100) {
			p1.dt = time - p1.t;
		}

		p1.x += p1.vx * p1.dt;
		p1.y += p1.vy * p1.dt;
		p1.z += p1.vz * p1.dt;
		p1.t += p1.dt;
		p1.dt = 0.0;
		time = p1.t;

		if ((p1.i_copy > -1) && (jm > -2)) {
			double delta = p1.t - particles[p1.i_copy].t;
			particles[p1.i_copy].x += particles[p1.i_copy].vx * delta;
			particles[p1.i_copy].y += particles[p1.i_copy].vy * delta;
			particles[p1.i_copy].z += particles[p1.i_copy].vz * delta;
			particles[p1.i_copy].t = p1.t;
			particles[p1.i_copy].dt = 0.0;
		}

		particles[im] = p1;

		need_virt_particle_retime = reform(im, jm);

		if (im >= NP) {
			if (particles[im].i_copy > -1)
				retime(im);
			if (jm > -2) {
				int e = im - NP;
				retime(e);
			}
		}
		else {
			retime(im);
			if ((particles[im].i_copy > -1) && ((need_virt_particle_retime == true) || (jm > -2) || ((jm < -10) && (jm != -100)))) {
				retime(particles[im].i_copy);
			}
		}

		if (jm >= 0) {
			++COLL_COUNT;

			if (jm >= NP) {
				if (particles[jm].i_copy > -1)
					retime(jm);
				int e = jm - NP;
				retime(e);
			}
			else {
				retime(jm);
				if (particles[jm].i_copy > -1)
					retime(particles[jm].i_copy);
			}
		}
	}

	particle *p = particles;
	for (i = 0; i < NP * 2; ++i, ++p)
		(*p).t -= time;
	Event *t = time_queue;
	++t;
	for (i = 1; i < last; ++i, ++t)
		(*t).t -= time;
	time = 0.0;
}


/*
������� ��������� ������� ��������� �������
����� ������� ������������� ������������ ������������� ���� ������
���������:
file_name - ��� ����� ��� ���������� ������
*/
void image(int steps, short accuracy, std::string file_name) {
	const int W = (int(L) + 1) * 2 * accuracy; // number of dots for all system
	int img[10000];
	double t_global, dt;

	printf("INFO: Image started for %d steps with accuracy %d\n", steps, accuracy);

	for (short g = 0; g < 10000; ++g) img[g] = 0;

	for (short h = 0; h < steps; ++h) {
		step();
		// ������� ������ � ����� ���������� � ���������
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("%d / %d", h + 1, steps);

		t_global = get_maximum_particle_time();

		for (int i = 0; i < NP; ++i) {
			dt = t_global - particles[i].t;
			int m = int((L + particles[i].x + particles[i].vx * dt) * 10.0 * accuracy) / 10;
			++img[m];
		}
	}

	double x = 0.0, delta_x = 1.0 / accuracy;
	FILE *save_file = fopen(file_name.c_str(), "w+");
	for (int f = 0; f < W; ++f) {
		fprintf(save_file, "%f %d\n", x, img[f]);
		x += delta_x;
	}
	fclose(save_file);

	printf("\n INFO: Image completed. Information saved to file: %s \n", file_name.c_str());
}


/*
������� �������� ������� ������� � ����
��������� ������������� ������������ ������ � ����,
������������ ��������� ������.
���������:
x1 - ��������� X ���������� "�����"
x2 - �������� X ���������� "�����"
file_name - ��� ����� ��� ���������� ������
*/
void profile(double x1, double x2, int dots_per_particle, int steps, std::string file_name) {
	double x, y, z, dt;
	int i, j;
	FILE *profile_file = fopen(file_name.c_str(), "w+");

	printf("INFO: Profile started.\n");

	for (i = 0; i < dots_per_particle; ++i) {
		for (j = 0; j < steps; ++j)
		    step();

		double t_global = get_maximum_particle_time();

		for (j = 0; j < NP; ++j) {
			dt = t_global - particles[j].t;
			x = L + particles[j].x + particles[j].vx * dt;
			if (x <= x2 && x >= x1) {
				y = A + particles[j].y + particles[j].vy * dt;
				z = A + particles[j].z + particles[j].vz * dt;
				// ��������� Y � Z ���������� �������
				fprintf(profile_file, "%.15le\n", y);
				fprintf(profile_file, "%.15le\n", z);
			}
		}
	}

	fclose(profile_file);

	printf("INFO: Profile completed. Information saved to file: %s\n", file_name.c_str());
}


/*
������� ������ �������, ��������� ������� ��� ��������� ������� �� �������� ���������
� �������� ������������ ����� �� ���������.
���������:
compress_to_etta - �������� ���������
delta_etta - ������������ ����������� ��� �� ���������
steps - ���������� ���������� �� ���� ������� � ������� �����
        ���������� �������� �� ���������
type - ��� ������:
       0 - ����������� ������ ����� ������
	   1 - ����������� ������ ������ ������
	   2 - ����������� ������������ ��� ������ �� ���������� ����������
*/
void compress(double compress_to_etta, double delta_etta, int steps, int type) {
	int m;
	double etta = (PI * NP) / (6.0 * A * A * (L - 1.0));
	double min = 1.0e+100, max = -1.0, x, dL, dx;
	double t_global, dt;
	double delta_L = fabs(L - 2.0 - ((PI*NP) / (etta + delta_etta)) / (6.0 * A * A) + 1.0);
	double L_ideal = ((PI * NP) / compress_to_etta) / (6.0 * A * A) + 1;

	printf("\n INFO: Start to change system density... \n");
	printf("etta = %.16le, should be equal to %.15le", etta, compress_to_etta);

	// ����� ��������� � ��������� � 12 ������
	while (fabs(etta - compress_to_etta) > 1.0e-12) {
		if (etta < compress_to_etta) {
			max = -100.0;
			min = 1.1e+10;

			// ������� ��� ������� � ������� ����� ������� ����� ��������������
			// ��� ������� � ����� ������� ���������� ������, �������� ������
			// ������������� �� ��������� ������, ����� ����� �� ������� �����
		    // ����������� ������ �� ���������� ������
			t_global = get_maximum_particle_time();

			for (int i = 0; i < NP; ++i) {
				dt = t_global - particles[i].t;
				x = L + particles[i].x + particles[i].vx * dt;
				if (x < min) min = x;
				if (x > max) max = x;

				if (particles[i].i_copy > 0) {
					m = particles[i].i_copy;
					dt = t_global - particles[m].t;
					x = L + particles[m].x + particles[m].vx * dt;
					if (x < min) min = x;
					if (x > max) max = x;
				}
			}

			min = min - 1.0;
			max = 2.0 * L - max - 1.0;

			dL = min;
			if (dL > max) dL = max;

			// ������� �� ������� � �������� � �� ������� ������
			dL = dL / 1.1;
			if (dL < 0.1e-8) dL = 0.01e-30;
			if (dL > delta_L) dL = delta_L;

			dx = dL / 2.0;

			// ���� ��������� �������� ������ ����� ������� ������ ��� ���������
			// �� ���������� ������ ������ �������� L, ����� ��������� ������� � ���������
			if (L - dL < L_ideal) {
				dL = 0.0;
				L = L_ideal;
			}
		}
		else {
			dL = L - L_ideal;
			if (dL < -delta_L) dL = -delta_L; // ��� ���������� �������
			dx = dL / 2.0;
		}

		if (type == 0) {
			// ������� ��� ������� �����, ����� ������� ���������� ������ ����� ������.
			// ���������� ������ �� ������ ������ �� ���������, �.�. ����� ����������� ������
			// ����� �� �������� ��� ������ �� �� �� ���������� - � ����� �� ������� ���
			// ������� ����� �� dL/2.0 � ���������� ��� ������ � ������ ������� �� dL/2.0.
			for (int w = 0; w < NP; ++w) {
				particles[w].x -= dx;
				if (particles[w].i_copy >= NP) particles[particles[w].i_copy].x -= dx;
			}
		}
		if (type == 1) {
			// ������� ������ ������ ������ (�� �� �� ����� ��� � ��� ����� ������, ��
			// � ���� ������ ��� ������� ��������� ������ �� dL/2.0)
			for (int w = 0; w < NP; ++w) {
				particles[w].x += dx;
				if (particles[w].i_copy >= NP) particles[particles[w].i_copy].x += dx;
			}
		}

		// �������� ������� - ���������� ���������� ��������� ��������� ���� ������
		L -= dL;

		// ��������� ��������� ������� � ����� ��������� ��� ���������� �����
		// ��������� � ������� ����������� �������������
		save("tmp");
		load_seed("tmp");

		for (short i = 0; i < steps; ++i)
			step();

		// ������������ ��������� ����� ������
		etta = (PI * NP) / (6.0 * A * A * (L - 1.0));

		// ������� ������ � ����� ���������� � ���������
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("etta = %.16le, should be equal to %.16le", etta, compress_to_etta);
	}
	printf("\n INFO: System density was sucessfully changed to %.16le \n", etta);
}


void init(std::string file_name) {
	using namespace std;
	clock_t start, end, result;
	char command[255], parameter[255];
	int i, steps;
	ifstream command_file(file_name.c_str());

	while (!command_file.eof()) {
		command_file.getline(command, 255, ' ');
		string str_command = command;

		printf("\n\n<==========================>\n");

		if (str_command.compare("new") == 0) {
			int NN;
			double etta;

			command_file >> NN;
			command_file >> etta;
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������
			new_seed(NN, etta);

			print_system_parameters();
		}
		if (str_command.compare("load") == 0) {
			command_file.getline(parameter, 255, '\n');
			load_seed(parameter);

			print_system_parameters();
		}
		if (str_command.compare("step") == 0) {
			command_file >> steps;
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������

			printf(" INFO: Step start. \n");

			start = clock();

			for (i = 0; i < steps; ++i) {
				step();
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
				printf("%d / %d", i + 1, steps);

				if (i % 50 == 0) {
					FILE *history_file = fopen("history.txt", "w+");
					fclose(history_file);
				}
			}

			end = clock();
			result = end - start;

			printf("\n INFO: finished %d collisions per particle \n", steps);
			printf("Total Time = %f seconds. \n", double(result / CLOCKS_PER_SEC));

			double Lx = 0.0;
			for (int i = 0; i < NP; i++) {
				Lx += (particles[i].y + A)*particles[i].vz - (particles[i].z + A)*particles[i].vy;
			}

			printf("\n Lx = %.15le\n", Lx);
		}
		if (str_command.compare("image") == 0) {
			command_file >> steps; // ����� ���������� �� ���� ������� �� ����� ���������
			command_file >> i; // ��������. ����� ����� ������� �� ���� ������ �������
			command_file.getline(parameter, 255, '\n');
			image(steps, i, parameter);
		}
		if (str_command.compare("profile") == 0) {
			int dots_for_each_particle;
			double x1, x2;
			command_file >> x1;
			command_file >> x2;
			command_file >> dots_for_each_particle;
			command_file >> steps;
			command_file.getline(parameter, 255, '\n');
			profile(x1, x2, dots_for_each_particle, steps, parameter);
		}
		if (str_command.compare("save") == 0) {
			command_file.getline(parameter, 255, '\n');
			save(parameter);
			printf("\n INFO: particles coordinates saved to '%s' \n", parameter);
		}
		if (str_command.compare("compress") == 0) {
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������
			printf("\n Sorry, 'compress' command was deprecated, we need to use \n");
			printf("'compress_two_walls'(short form: 'compresst') instead. You can also use ");
			printf("'compress_left_wall'(short form: 'compressl') \n or 'compress_right_wall'");
			printf("(short form: 'compressr') commands to change system density.\n");
			exit(1);
		}
		if ((str_command.compare("compress_two_walls") == 0) || (str_command.compare("compresst") == 0)) {
			double etta, delta_etta;
			command_file >> etta; // ��������� ���������
			command_file >> delta_etta;  // ����������� ���������� �������� ��������� ���������
			command_file >> steps;  // ���������� ���������� ����� ������� ���� ������
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������
			compress(etta, delta_etta, steps, 2);

			print_system_parameters();
		}
		if ((str_command.compare("compress_left_wall") == 0) || (str_command.compare("compressl") == 0)) {
			double etta, delta_etta;
			command_file >> etta; // ��������� ���������
			command_file >> delta_etta;  // ����������� ���������� �������� ��������� ���������
			command_file >> steps;  // ���������� ���������� ����� ������� ���� ������
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������
			compress(etta, delta_etta, steps, 0);

			print_system_parameters();
		}
		if ((str_command.compare("compress_right_wall") == 0) || (str_command.compare("compressr") == 0)) {
			double etta, delta_etta;
			command_file >> etta; // ��������� ���������
			command_file >> delta_etta;  // ����������� ���������� �������� ��������� ���������
			command_file >> steps;  // ���������� ���������� ����� ������� ���� ������
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������
			compress(etta, delta_etta, steps, 1);

			print_system_parameters();
		}
		if (str_command.empty()) break;
	}
	FILE *result_flag = fopen("result", "w+");
	fclose(result_flag);
}


int main()
{
	FILE *history_file = fopen("history.txt", "w+");
	fclose(history_file);


	init("program.txt");
	return 0;
}