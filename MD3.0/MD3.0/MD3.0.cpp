/*
������������� ������
*/
#include "stdafx.h"  // ����������� ����������
#include <time.h>    // ���������� ��� ��������� �������
#include <fstream>   // ���������� ��� ������ � �������
#include <cmath>     // ���������� ��� ������ �������������� �������


/*
���������� ���������
*/

// ���������� ����� �� Y, Z � X(K2), �� ������� ����������� �����.
// ���������� ����� ������������� ����������� ��� �������� �������
// ��� ��� ������ (��. ������� load).
short K, K2;

// ����� ������ �� ��� ������, ���������������� � �������� load � new_seed
int NP = 6976;

// ���������� ������� ������������ � �������
int COLL_COUNT = 0;

#define PI 3.141592653589793238462

// ��������� ������, �������� � load()
double A, dA, L, dL;

// ���������� ���������� ��� �������� ����� ������������ ������� ���� ������
double global_E = 0.0;

// ������ ���������� �������� � ������� �������.
int last;

// ������ "�������"
typedef struct Event_ {
	double t;
	int im, jm, vp1, vp2;
} Event;

// ������� ������� - ���������� 8192 ��������
// (��� ������ ���� �����-������� ������, ������� ��� ������������ ����� ������)
Event time_queue[8192];

// ������ "�������"
// x, y, z - ���������� �������
// vx, vy, vz - �������� ��������� �������
// t - ����������� ����� �������
// dt - ����� �� ���������� ������� ���� �������
// x_box, y_box, z_box - ����� ������, � ������� ��������� �������
// ti - ����� ������� ������� � ������ �������
// box_i - ����� ������� � ������
typedef struct particle_ {
	double x[4], y[4], z[4], vx, vy, vz, t, dt;
	int x_box[4], y_box[4], z_box[4], ti, box_i[4];
} particle;

// ������ ������, ������ ������� N + ����������
// � ������� ������� � ����� ������ ������� ������.
particle particles[8192];

// ������. ����� ������� ������� �� ��������� ������,
// ������ ������ �������� � ���� ��������� ����������� ������
// x1, y1, z1, x2, y2, z2 - ���������� ����� � ������ ������ ������
// particles[] - ������ ���� ������, ����������� � ������ ������
// end - ������ ��������� ������� � ������ ������ ������ ������
typedef struct Box_ {
	double x1, y1, z1, x2, y2, z2;
	int particles[32];
	short end;
} Box;

// ������ ������ ��� ����� ������
Box boxes_yz[16][16][64];


// ������ � �������� ������, ��� ������� ���� ��������� ������� �������
// ������������ �� ������ ������� ��������� ��� ���������� ������� �������
// ��������� ������
int particles_for_check[100];
int particles_for_check_count = 0;

/*
��� ������� ������� �� ����� ��������� ������� �������:
A - ������ ������� (������� ������, ��� ����� ���������� ����� �������) �� Y � Z.
L
N - ����� ������ � �������
etta - ������� ������������� ��������� �������
*/
void print_system_parameters() {
	long double etta = (4.0 * PI * NP) / (3.0 * A * A * (L - 2.0));
	printf("\n\n| A    | %.15le |\n| L    | %.15le |\n", A, L);
	printf("| N    | %d                  |\n", NP);
	printf("| etta | %.15le |\n", etta);
}

/*
������� ���������� ���������� ����������� ����� ������ � �������,
��� ��������� �������������� ��� ������� �� �������
*/
double get_maximum_particle_time() {
	double t_max = -1.0e+20;

	for (int i = 0; i < NP; ++i) {
		if (particles[i].t > t_max) {
			t_max = particles[i].t;
		}
	}

	return t_max;
}

double get_minimum_particle_time() {
	double t_min = 1.0e+20;

	for (int i = 0; i < NP; ++i) {
		if (particles[i].t < t_min)
			t_min = particles[i].t;
	}

	return t_min;
}


/*
������� ��� �������� ��������� �������, � ��� �� ���������, ���
��� ������� ��������� ������ �������, ��� ������ ������� � ��� ����������
��������� �������, ������� ��������� � ���������� ������� ������� � ��.

� ������ ������������� ����� ������� ������ ������� ���������� ������ ���������
� ������� �������������� ���������� �� ������������ ��������.

������� �������� ��������, ���������� ������������ � ����� ��������
��������� � ���������.
*/
int check_particles() {
	double E = 0.0;
	double dt = 0.0;
	double t_global = get_maximum_particle_time();

	for (int i = 0; i < K + 1; i++) {
		for (int j = 0; j < K + 1; j++) {
			for (int f = 0; f < K2 + 1; f++) {
				if (boxes_yz[i][j][f].end > 30) {
					printf("!! %d %d %d  end = %d", i, j, f, boxes_yz[i][j][f].end);
					throw "AAA";
				}

				for (int p = 0; p < boxes_yz[i][j][f].end + 1; p++) {
					if (boxes_yz[i][j][f].particles[p] > NP * 4) {
						printf("\n !! %d %d %d  boxes_yz[i][j][f].particles[p] %d \n", i, j, f, boxes_yz[i][j][f].particles[p]);
						throw "AAA";
					}
				}
			}
		}
	}

	for (int i = 0; i < NP; i++) {
		particle p1 = particles[i];
		Box p1_box = boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]];

		dt = t_global - p1.t;
		p1.x[0] += p1.vx * dt;
		p1.y[0] += p1.vy * dt;
		p1.z[0] += p1.vz * dt;

		E += p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;

		for (int j = i + 1; j < NP; j++) {
			particle p2 = particles[j];

			dt = t_global - p2.t;
			p2.x[0] += p2.vx * dt;
			p2.y[0] += p2.vy * dt;
			p2.z[0] += p2.vz * dt;

			double dx = p1.x - p2.x;
			double dy = p1.y - p2.y;
			double dz = p1.z - p2.z;
			double r = dx*dx + dy*dy + dz*dz;

			if ((r < 4.0) && (4.0 - 1.0e-5 - r > 0.0)) {
				printf("\n\n particles %d %d r = %.15le", i, j, r);
				throw "Particles overlaps!";
			}
		}

		if (p1.x[1] < 0.0 && (p1.y[0] < 2.0 || p1.y[0] > A - 2.0)) {
			printf("\n\n %d p1.y[0] %.15le", i, p1.y[0]);
			throw "aa";
		}
		if (p1.x[2] < 0.0 && (p1.z[0] < 2.0 || p1.z[0] > A - 2.0)) {
			printf("\n\n %d p1.z[0] %.15le", i, p1.z[0]);
			throw "aa";
		}

		// ��������� ������� ����� ��� ���� ������
		if ((p1.x_box[0] > K2) || (p1.y_box[0] > K) || (p1.z_box[0] > K) ||
			(p1.x_box[0] < 0) || (p1.y_box[0] < 0) || (p1.z_box[0] < 0)) {
			throw "Particle locates in incorrect cell.";
		}

		// ��������� ��� ������� ��������� � ������
		if ((p1.x[0] > L && p1.x[0] < 0.0) ||
			(p1.y[0] > A && p1.y[0] < 0.0) ||
			(p1.z[0] > A && p1.z[0] < 0.0)) {
			printf(" \n Particle %d, %.15le, %.15le, %.15le \n ", i, p1.x[0], p1.y[0], p1.z[0]);
			throw "Particle is out of the system boundaries.";
		}

		// ��������� ��� ������� ��������� � ���������� �������
		if (((p1.x[0] < p1_box.x1) && (abs(p1.x[0] - p1_box.x1) > 1.0e-14)) ||
			((p1.x[0] > p1_box.x2) && (abs(p1.x[0] - p1_box.x2) > 1.0e-14)) ||
			((p1.y[0] < p1_box.y1) && (abs(p1.y[0] - p1_box.y1) > 1.0e-14)) ||
			((p1.y[0] > p1_box.y2) && (abs(p1.y[0] - p1_box.y2) > 1.0e-14)) ||
			((p1.z[0] < p1_box.z1) && (abs(p1.z[0] - p1_box.z1) > 1.0e-14)) ||
			((p1.z[0] > p1_box.z2) && (abs(p1.z[0] - p1_box.z2) > 1.0e-14))) {

			printf("Vilet za granicy %d \n", i);

			printf("Granizy:\n");
			printf("X : [%.15le ; %.15le]\n", p1_box.x1, p1_box.x2);
			printf("Y : [%.15le ; %.15le]\n", p1_box.y1, p1_box.y2);
			printf("Z : [%.15le ; %.15le]\n", p1_box.z1, p1_box.z2);
			printf("x, y, z: %.15le %.15le %.15le\n", p1.x[0], p1.y[0], p1.z[0]);

			printf("p1.t = %.15le, p.im = %d, p1.jm = %d \n", p1.t, time_queue[p1.ti].im, time_queue[p1.ti].jm);

			throw "Particle is out of the cell boundary.";
		}

		
		// ��������� ��� ������� �������� � ����� �� ����� � �������
		bool w = false;
		for (int ty = 0; ty <= boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].end; ++ty) {
			if (boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].particles[ty] == i) w = true;
		}
		if (w == false) {
			for (int t = 0; t <= boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].end; ++t) {
				printf(" %d ", boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].particles[t]);
			}

			printf("\n Particle %d doesn't store in the cell.", i);

			throw "Particle doesn't store in the cell.";
		}

		/*
		// ��������� �� ���������� �� ������� � ������� ����� ��������� �������
		if (time_queue[p1.ti].im != i && time_queue[p1.ti].jm != i) {
		Event e = time_queue[p1.ti];
		printf("\n i = %d ; im = %d ; jm = %d ; ti = %d ", i, e.im, e.jm, p1.ti);
		throw "Particle has no correct link to the event.";
		}
		*/
	}

	// ��������� ������� �������� ���������� ������������ ������� ������� ��
	// ��������� �������, ������� ���� ��� �������� ������� �� ����� � load
	if (abs(E - global_E) > global_E*0.1e-8) {
		printf("\nENERGY was changed: \n E_seed = %.15le \n E_now= %.15le \n", global_E, E);
		throw "ENERGY was changed.";
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
	int j = i >> 1;  // ��� ����� ������, �� �� ����� ��� j = i/2, ������ �������
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
i - �������, � �������� ������� ��������� ����� �������
j - ����� ������� � ������� ��������� ������� i ��� ����� �������
���������� �� ������� ��� ����������� ����� �������������
��������� ������� ��� ����� �������� �������
*/
void add_event(int &i, int &j, int &vp1, int &vp2) {
	/*
	������������ ������ ����� ������ ������� �� ������ �������
	����������� ������� �������, ����� ������� �� �������� ����� t,
	��������� ������� �� ����� ���������� ����� �� ������� � �������
	��������� ������
	*/
	double t = particles[i].t + particles[i].dt;

	/*
	������� ������� � ������ ����� ��� ������ �������.
	���������� �������� ��� ������� ���� ������ � ���������
	��� ��������� �� ������, ���� ������ ������� ��������� ������ ���
	������ �������
	*/
	particles[i].ti = get_up(last, t);

	/*
	���� ����� ������� - ��� ������� ������������ ���� ������, ��
	���������� ��� ������ ������� ��������� ������ � � ����� �������
	*/
	if (j >= 0) {
		particles[j].dt = particles[i].dt;
		particles[j].ti = particles[i].ti;
	}

	// ���������� ����� ������� � ��������� ������ � ������ �����
	time_queue[particles[i].ti].im = i;
	time_queue[particles[i].ti].jm = j;
	time_queue[particles[i].ti].vp1 = vp1;
	time_queue[particles[i].ti].vp2 = vp2;
	time_queue[particles[i].ti].t = t;

	// ����������� ����� ������� �� 1
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


/*
������� �������� ������� ������� �� ������� �������

���������:
i - ����� ������� ��� ������, ������� �������� ���������� �������
*/
void clear_particle_events(int &i) {
	// particles[i].ti - ������ �� ������ �������, ������� ������ ���� �������
	int f = particles[i].ti;

	if (f > 0) {
		int e = -100, pv1 = 0, pv2 = -1;
		int kim = time_queue[f].im;
		int kjm = time_queue[f].jm;

		if (time_queue[f].im == i) {
			/*
			���� �� ������� ������� ������������ ���� ������, ��
			���������� ����������� ������ ������� � �� �� �����,
			� ������� ��������� �������, ��� ������� �� �������
			��� ������� � �������� ������� ������������ ���� ������
			�� ������� "-100", ����� ������ ������� ������ ������� ��
			����� ��������������� ������������ � ����� ����� ��� ��
			����� ���������� ����� �������.
			*/
			if (time_queue[f].jm >= 0) {
				double dt = particles[kim].t - particles[kjm].t;  // ������� � ������� �������

				for (int r = 0; r < 4; ++r) {
					if (particles[kjm].x[r] > 0.0) {
						particles[kjm].x[r] += particles[kjm].vx*dt;
						particles[kjm].y[r] += particles[kjm].vy*dt;
						particles[kjm].z[r] += particles[kjm].vz*dt;
					}
				}
				particles[kjm].t = particles[kim].t;
				particles[kjm].dt = (particles[kjm].dt - dt) / 1.1;

				delete_event(f);
				add_event(kjm, e, pv1, pv2);
			}
			else delete_event(f);
		}
		else if (time_queue[f].jm == i) {
			double dt = particles[kjm].t - particles[kim].t;  // ������� � ������� �������

			for (int r = 0; r < 4; ++r) {
				if (particles[kim].x[r] > 0.0) {
					particles[kim].x[r] += particles[kim].vx*dt;
					particles[kim].y[r] += particles[kim].vy*dt;
					particles[kim].z[r] += particles[kim].vz*dt;
				}
			}
			particles[kim].t = particles[kjm].t;
			particles[kim].dt = (particles[kim].dt - dt) / 1.1;

			delete_event(f);
			add_event(kim, e, pv1, pv2);
		}
		particles[i].ti = 0;
	}
}


/*
������� ������� ���������� ������� ��� �������

���������:
i - ����� �������, ��� ������� �� ������ ���������� ��������� �������
*/
void retime(int &i) {
	particle p1 = particles[i];
	Box p1_box = boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]];
	int jm = -1;  // ���������� ��� ���������� ���� ���������� �������
	double dt, dt_min = 1.0e+20;  // ���������� ��� �������� ������� ��� ���������� �������

	if (p1.vx < 0.0) {
		dt_min = (p1_box.x1 - p1.x[0]) / p1.vx;
		jm = -2;  // ������� ����������� ������� �1 ������, � ������� ��������� ������� 

	    // ���� �� ��������� ������ ��������� ������, �� ���������� ����� ����������
		// � ��������� �������
		if (p1.x_box[0] == 1) {
			dt_min = (p1_box.x1 + 1.0 - p1.x[0]) / p1.vx;
			jm = -1;  // ������� ������������ � ��������� �������
		}
	}
	else {
		dt_min = (p1_box.x2 - p1.x[0]) / p1.vx;
		jm = -4;  // ������� ����������� ������� �2 ������, � ������� ��������� �������

		// ���� �� ��������� ������ ��������� ������, �� ���������� ����� ����������
		// � ��������� �������
		if (p1.x_box[0] == K2 - 1) {
			dt_min = (p1_box.x2 - 1.0 - p1.x[0]) / p1.vx;
			jm = -1;  // ������� ������������ � ��������� �������
		}
	}

	if (p1.vy < 0.0) {
		dt = (p1_box.y1 - p1.y[0]) / p1.vy;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -5;  // ������� ����������� ������� Y1 ������, � ������� ��������� �������  
		}
		if (p1.y_box[0] == 1) {
			dt = (p1_box.y1 + 2.0 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -12;  // ������� �������� ������ �������
			}
		}
		if (p1.y_box[0] == K-1) {
			dt = (p1_box.y2 - 2.0 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -11;  // ������� �������� ������ �������
			}
		}
	}
	else {
		dt = (p1_box.y2 - p1.y[0]) / p1.vy;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -6;  // ������� ����������� ������� Y2 ������, � ������� ��������� �������
		}
		if (p1.y_box[0] == K - 1) {
			dt = (p1_box.y2 - 2.0 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -11;  // ������� �������� ������ �������
			}
		}
		if (p1.y_box[0] == 1) {
			dt = (p1_box.y1 + 2.0 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -12;  // ������� �������� ������ �������
			}
		}
	}

	if (p1.vz < 0.0) {
		dt = (p1_box.z1 - p1.z[0]) / p1.vz;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -7;  // ������� ����������� ������� Z1 ������, � ������� ��������� �������
		}
		if (p1.z_box[0] == 1) {
			dt = (p1_box.z1 + 2.0 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -14;  // ������� �������� ������ �������
			}
		}
		if (p1.z_box[0] == K - 1) {
			dt = (p1_box.z2 - 2.0 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -13;  // ������� �������� ������ �������
			}
		}
	}
	else {
		dt = (p1_box.z2 - p1.z[0]) / p1.vz;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -8;  // ������� ����������� ������� Z2 ������, � ������� ��������� �������
		}
		if (p1.z_box[0] == K - 1) {
			dt = (p1_box.z2 - 2.0 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -13;  // ������� �������� ������ �������
			}
		}
		if (p1.z_box[0] == 1) {
			dt = (p1_box.z1 + 2.0 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -14;  // ������� �������� ������ �������
			}
		}
	}

	double temp, dx, dy, dz, dvx, dvy, dvz, d, dv, bij, dr;
	int s, n, r, q, w, vp1 = 0, vp2 = -1;

	for (int iv = 0; iv < 4; iv++) {
		if (p1.x[iv] > 0.0) {
			// �������� �� �������, ��������� � ������, � ������� ��������� ������� i
			for (r = p1.x_box[iv] - 1; r < p1.x_box[iv] + 2; r++)
				for (q = p1.y_box[iv] - 1; q < p1.y_box[iv] + 2; q++)
					for (w = p1.z_box[iv] - 1; w < p1.z_box[iv] + 2; w++) {

						// ���� ������ ������ ������� �� ������� ������� ��
						// ��������� �� ��������� ��� �����
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

						// �������� �� ���� �������� � ��������� ������ � ��������� �����������
						// ������������ ���� ������ � �������� i
						for (s = 0; s <= boxes_yz[q][w][r].end; ++s) {
							n = boxes_yz[q][w][r].particles[s];
							int jv = 0;

							if (n >= NP) {
								jv = n / NP;
								n = n % NP;
							}

							/*
							C�������� � ���������� p ��� ������ ������� n,
							����� ����� ���������� ��� �������� ������ � ��� ���������
							� ���������� ������
							*/
							particle p = particles[n];

							/*
							������������ ������� ������������ ������� ���� ������,
							��� ���������� ����� �������������� �� ����� ����� � ������������
							����������� � ����� ������������ � �������, ��� ��� ������� �����
							����� ���������� ����������� �����
							*/
							temp = p1.t - p.t;

							dvx = p.vx - p1.vx;
							dvy = p.vy - p1.vy;
							dvz = p.vz - p1.vz;

							dx = p.x[jv] + p.vx * temp - p1.x[iv];
							dy = p.y[jv] + p.vy * temp - p1.y[iv];
							dz = p.z[jv] + p.vz * temp - p1.z[iv];

							bij = dx * dvx + dy * dvy + dz*dvz;

							if (bij < 0.0) {
								dv = dvx * dvx + dvy * dvy + dvz*dvz;
								dr = 4.0 - dx * dx - dy * dy - dz * dz;
								// ������������ ������������ � ��������� ��� ���������� ������� ����������
								d = bij * bij + dv * dr;

								// ���� ������������ ������ ���� �� ���������� ��������
								if (d > 0.0) {
									dt = -(sqrt(d) + bij) / dv;

									/*
									� ������� � ����������� ������� ������ ���������� ����� �� �� ����������,
									� ���������� �������� ����� dt ����� ������� ��� ���������� ��������� ���
									������� p, ����� ������� �� �������� ����������� �������� ����� ������ ���
									������������� ���������� �� �������� ���������� ������� ��� ������� p �
									������� ������� - ������ �� �� ������������ ������� ��� ������� p ��� ���.
									*/
									temp += dt;

									if ((dt < dt_min) && (dt >= 0.0) &&
										((temp < p.dt) || ((abs(temp - p.dt) < 0.1E-15) &&
											(time_queue[p.ti].im == n) &&
											(time_queue[p.ti].jm == -100)))) {
										dt_min = dt;
										jm = n;
										vp1 = iv;
										vp2 = jv;
									}
								}
							}
						}
					}
		}
	}

	/*
	���� ����� ������� - ��� ������������ ���� ������, �� ��� ������ �������
	���������� ������������ � ���������� ��������� ������� � ��������������
	��� ������� �� ������� � �������� i
	*/
	if (jm >= 0) {
		dt = p1.t - particles[jm].t;

		if (dt < 0.0 && dt > -1.0e-15)
			dt = 0.0;

		particles[jm].t = p1.t;
		particles[jm].dt = dt_min;

		for (int j = 0; j < 4; j++) {
			if (particles[jm].x[j] > 0.0) {
				particles[jm].x[j] += particles[jm].vx * dt;
				particles[jm].y[j] += particles[jm].vy * dt;
				particles[jm].z[j] += particles[jm].vz * dt;
			}
		}

		clear_particle_events(jm);
	}

	particles[i].dt = dt_min;

	add_event(i, jm, vp1, vp2);

	/*
	for (int h = 0; h < particles_for_check_count; h++) {
		if (i == particles_for_check[h] || jm == particles_for_check[h] || i == NP + particles_for_check[h] || jm == NP + particles_for_check[h]) {
			FILE *history_file = fopen("history.txt", "a");
			fprintf(history_file, "\n\n retime result %d:  %d, dt = %.15le\n", i, jm, particles[i].dt);
			fprintf(history_file, "particle %d p.x = %.15le, p.y = %.15le, p.z = %.15le\n", i, particles[i].x, particles[i].y, particles[i].z);
			fclose(history_file);
		}
	}
	*/

	/*
	� ������ ���� ��� ������� ���������� ������� �� �������� ������������� �����
	���������� �������� ���������� ��������� � ����������� ���������� ����������
	� ������������ �������
	*/
	if (dt_min < -1.0e-11) {
		printf("\n retime result: %d %d, %.16le\n ", i, jm, dt_min);
		printf("\n p1.x = %.16le, p1.y = %.16le, p1.z = %.16le \n", p1.x[0], p1.y[0], p1.z[0]);
		printf("\n p1.vx = %.16le, p1.vy = %.16le, p1.vz = %.16le \n", p1.vx, p1.vy, p1.vz);
		printf("\n p1.x_box = %d, p1.y_box = %d, p1.z_box = %d", p1.x_box[0], p1.z_box[0], p1.y_box[0]);
		printf("\n im = %d, jm = %d, dt = %.16le, A = %.16le", i, jm, dt_min, A);
		printf("\n p1.box.x = [%.16le; %.16le]", boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].x1, boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].x2);
		printf("\n p1.box.y = [%.16le; %.16le]", boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].y1, boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].y2);
		printf("\n p1.box.z = [%.16le; %.16le]", boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].z1, boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].z2);
		throw "dt < 0!";
	}
}


void change_particles(int i, int vp1, int vp2) {
	Box box1, box2;
	int y_box1, z_box1, x_box1, y_box2, z_box2;
	int j;

	x_box1 = particles[i].x_box[vp1];
	y_box1 = particles[i].y_box[vp1];
	z_box1 = particles[i].z_box[vp1];

	box1 = boxes_yz[y_box1][z_box1][x_box1];
	box1.particles[particles[i].box_i[vp1]] = i + NP * vp2;
	boxes_yz[y_box1][z_box1][x_box1] = box1;

	j = particles[i].box_i[vp1];
	particles[i].box_i[vp1] = particles[i].box_i[vp2];

	y_box2 = particles[i].y_box[vp2];
	z_box2 = particles[i].z_box[vp2];

	box2 = boxes_yz[y_box2][z_box2][x_box1];
	box2.particles[particles[i].box_i[vp2]] = i + NP * vp1;
	boxes_yz[y_box2][z_box2][x_box1] = box2;

	particles[i].box_i[vp2] = j;

	particles[i].y_box[vp1] = y_box2;
	particles[i].z_box[vp1] = z_box2;
	particles[i].y_box[vp2] = y_box1;
	particles[i].z_box[vp2] = z_box1;
}

/*
������� ������ ������� � "������" ��� ����������� �������� �������������
��������� �������.

���������:
im - ����� �������, ������������ ������������� ��������� �������
jm - ����� �������, ����� ������� �������� ����� �������
*/
void change_with_virt_particles(int &im, int &jm) {
	double k;

	if (jm > -7) {

		k = particles[im].y[0];
		particles[im].y[0] = particles[im].y[1];
		particles[im].y[1] = k;

		change_particles(im, 0, 1);

		if (particles[im].x[3] > 0.0) {
			k = particles[im].y[2];
			particles[im].y[2] = particles[im].y[3];
			particles[im].y[3] = k;

			change_particles(im, 2, 3);
		}
	}
	else {
		k = particles[im].z[0];
		particles[im].z[0] = particles[im].z[2];
		particles[im].z[2] = k;

		change_particles(im, 0, 2);

		if (particles[im].x[1] > 0.0) {
			k = particles[im].y[1];
			particles[im].y[1] = particles[im].y[3];
			particles[im].y[3] = k;

			change_particles(im, 1, 3);
		}
	}

}


/*
������� �������� ������ "������".

���������:
i - ����� �������, ����� ������� ���������� �������
need_to_check - ����, ����������� �� ��������� ����� �� ����� ��� ������
������� ��� ���, � ����� ��������� ����� (����� �� �������,
��� ����� ���������� �������).
*/
void create_virt_particle(int &im, int &jm) {
	Box box1;

	switch (jm)
	{
		case -11:     // ������� �� ������� Y = �-2
			particles[im].y[1] = -2.0;
			particles[im].x[1] = particles[im].x[0];
			particles[im].z[1] = particles[im].z[0];
			particles[im].x_box[1] = particles[im].x_box[0];
			particles[im].y_box[1] = 0;
			particles[im].z_box[1] = particles[im].z_box[0];

			box1 = boxes_yz[0][particles[im].z_box[1]][particles[im].x_box[1]];
			box1.end++;
			box1.particles[box1.end] = im + NP;
			boxes_yz[0][particles[im].z_box[1]][particles[im].x_box[1]] = box1;
			particles[im].box_i[1] = box1.end;

			break;
		case -12:     // ������� �� ������� Y = 2
			particles[im].y[1] = A + 2.0;
			particles[im].x[1] = particles[im].x[0];
			particles[im].z[1] = particles[im].z[0];
			particles[im].x_box[1] = particles[im].x_box[0];
			particles[im].y_box[1] = K;
			particles[im].z_box[1] = particles[im].z_box[0];

			box1 = boxes_yz[K][particles[im].z_box[1]][particles[im].x_box[1]];
			box1.end++;
			box1.particles[box1.end] = im + NP;
			boxes_yz[K][particles[im].z_box[1]][particles[im].x_box[1]] = box1;
			particles[im].box_i[1] = box1.end;

			break;
		case -13:     // ������� �� ������� Z = A-2
			//if (im == 619) printf("\n &&&&&&& \n");

			particles[im].z[2] = -2.0;
			particles[im].x[2] = particles[im].x[0];
			particles[im].y[2] = particles[im].y[0];
			particles[im].x_box[2] = particles[im].x_box[0];
			particles[im].z_box[2] = 0;
			particles[im].y_box[2] = particles[im].y_box[0];

			box1 = boxes_yz[particles[im].y_box[2]][0][particles[im].x_box[2]];
			box1.end++;
			box1.particles[box1.end] = im + NP*2;
			boxes_yz[particles[im].y_box[2]][0][particles[im].x_box[2]] = box1;
			particles[im].box_i[2] = box1.end;

			break;
		case -14:     // ������� �� ������� Z = 2
			particles[im].z[2] = A + 2.0;
			particles[im].x[2] = particles[im].x[0];
			particles[im].y[2] = particles[im].y[0];
			particles[im].x_box[2] = particles[im].x_box[0];
			particles[im].z_box[2] = K;
			particles[im].y_box[2] = particles[im].y_box[0];

			box1 = boxes_yz[particles[im].y_box[2]][K][particles[im].x_box[2]];
			box1.end++;
			box1.particles[box1.end] = im + NP*2;
			boxes_yz[particles[im].y_box[2]][K][particles[im].x_box[2]] = box1;
			particles[im].box_i[2] = box1.end;

			break;
	}

	// ���� � ������� ���� ��� ������, �� ���������� ������� ������ �����
	if ((particles[im].x[1] > 0) && (particles[im].x[2] > 0))
	{
		particles[im].x[3] = particles[im].x[0];
		particles[im].y[3] = particles[im].y[1];
		particles[im].z[3] = particles[im].z[2];
		particles[im].x_box[3] = particles[im].x_box[0];
		particles[im].z_box[3] = particles[im].z_box[2];
		particles[im].y_box[3] = particles[im].y_box[1];

		int y = particles[im].y_box[3];
		int z = particles[im].z_box[3];
		box1 = boxes_yz[y][z][particles[im].x_box[3]];
		box1.end++;
		box1.particles[box1.end] = im + NP*3;
		boxes_yz[y][z][particles[im].x_box[3]] = box1;
		particles[im].box_i[3] = box1.end;
	}
}


/*
������� �������� ���������� � ��������� ������� �� �����.
��� �������� ������� �� ��������������� ��������� ����������
� ������ ������� � ������� ���� �������.

���������:
loading_file - ������ �� ����, �� �������� ���������� �������
������ � ������� ��� ������.
*/
void load_information_about_one_particle(FILE *loading_file) {
	double a1, a2, a3;
	int i, x_box, y_box, z_box, end;

	// �������� ����� �������
	fscanf(loading_file, "%d\n", &i);

	// ��������� ���������� �������
	fscanf(loading_file, "%le %le %le\n", &a1, &a2, &a3);
	particles[i].x[0] = a1;
	particles[i].y[0] = a2;
	particles[i].z[0] = a3;

	// ��������� ���������� � ��������� �������
	fscanf(loading_file, "%le %le %le\n", &a1, &a2, &a3);
	particles[i].vx = a1;
	particles[i].vy = a2;
	particles[i].vz = a3;

	particles[i].t = 0.0;
	particles[i].dt = 0.0;
	particles[i].ti = 1;

	if (particles[i].y[0] < 2.0) {
		particles[i].x[1] = particles[i].x[0];
		particles[i].y[1] = particles[i].y[0] + A;
		particles[i].z[1] = particles[i].z[0];
	}
	if (particles[i].y[0] > A - 2.0) {
		particles[i].x[1] = particles[i].x[0];
		particles[i].y[1] = particles[i].y[0] - A;
		particles[i].z[1] = particles[i].z[0];
	}
	if (particles[i].z[0] < 2.0) {
		particles[i].x[2] = particles[i].x[0];
		particles[i].z[2] = particles[i].z[0] + A;
		particles[i].y[2] = particles[i].y[0];
	}
	if (particles[i].z[0] > A - 2.0) {
		particles[i].x[2] = particles[i].x[0];
		particles[i].z[2] = particles[i].z[0] - A;
		particles[i].y[2] = particles[i].y[0];
	}
	if (particles[i].x[1] > 0.0 && particles[i].x[2] > 0.0) {
		particles[i].x[3] = particles[i].x[0];
		particles[i].z[3] = particles[i].z[2];
		particles[i].y[3] = particles[i].y[1];
	}

	for (int j = 0; j < 4; j++) {
		if (particles[i].x[j] > 0.0) {
			// ���������� � ����� ������ ���������� �������� ����� �������
			x_box = short((particles[i].x[j] + dL) / dL);
			y_box = short((particles[i].y[j] + dA) / dA);
			z_box = short((particles[i].z[j] + dA) / dA);

			if (y_box < 0) y_box = 0;
			if (z_box < 0) z_box = 0;
			if (y_box > K) y_box = K;
			if (z_box > K) z_box = K;

			particles[i].x_box[j] = x_box;
			particles[i].y_box[j] = y_box;
			particles[i].z_box[j] = z_box;

			// ���������� ������� � ������
			end = particles[i].box_i[j] = ++boxes_yz[y_box][z_box][x_box].end;
			boxes_yz[y_box][z_box][x_box].particles[end] = i + NP*j;

			/*
			���� � ����� ������ ��������� ������� ����� ������,
			�� ��������� ���������, ������ ��������� �� ������.
			*/
			if (end > 10) {
				printf("Error: Too many particles in one box");
				printf("\n x_box, y_box, z_box = %d %d %d , end = %d \n", x_box, y_box, z_box, end);
				throw "Error: Too many particles in one box";
			}
		}
	}

	global_E += a1*a1 + a2*a2 + a3*a3;

	/*
	��������� ��� � ������ ������ ��� ������������� ������� ������
	*/
 	Box b = boxes_yz[particles[i].y_box[0]][particles[i].z_box[0]][particles[i].x_box[0]];
	for (short t = 0; t < particles[i].box_i[0]; ++t)
		for (short d = t + 1; d <= particles[i].box_i[0]; ++d) {
			if (b.particles[t] == b.particles[d])
			{
				printf("\n Duplicated particles in one box: %d \n", i);

				for (int j = 0; j < b.end; ++j)
					printf("%d ", b.particles[j]);

				throw "Duplicated particles in one box";
			}
		}
}


/*
������� �������� ������ � ������� �� �����.
����������� ��������� ������, ���������� � �������� ���� ������.

���������:
file_name - ��� �����, ����������� ��� ������ � �������.
*/
void load_seed(std::string file_name) {
	double a1, x, y, z;
	int i;
	FILE *loading_file;

	global_E = 0.0;

	/*
	��������� ��������� ���� ��� ������,
	����� �� ����� ���������� �� ������ � ������� �� ��������� ����
	���������� ����� �� ���������� ��� ��� ���������� ������� ��� ������.
	*/
	try {
		loading_file = fopen(file_name.c_str(), "r");
	}
	catch (...) {
		printf("\n ERROR: Can't read file '%s' - is it exist? \n", file_name.c_str());
		exit(1);
	}

	// ��������� ���������� ������ � �������
	fscanf(loading_file, "%i\n", &i);
	NP = i;
	// ��������� �������� A
	fscanf(loading_file, "%le\n", &a1);
	A = a1;
	// ��������� �������� L 
	fscanf(loading_file, "%le\n", &a1);
	L = a1;

	// EN: search for the appropriate values for dA, dL, K, K2
	// RU: ���� ���������� �������� ��� dA, dL, K, K2
	dA = 2.2;
	dL = 2.2;
	K = short(A / dA) + 1;  // ���������� ����� �� Y � Z
	K2 = short(L / dL) + 1; // ���������� ����� �� X
	dA = A / (K - 1);       // ������ ������� ����� �� X, Y � Z
	dL = L / (K2 - 1);      // ��������� ���, �����, ������ ���������
							// � ��������� ������

	/*
	   ����� ���� ��� �� ���������� ���������� ����� � �� �������
	   ���������� ���������������� ������ ������ ������, �������
	   ������� ������ ������ �� X, Y, Z, �.�. ��� ������ �����
	   �������������� ��� �������� ������� �� ����������� ������
	   ����� ���������, ������� ��������� � ���� �������.
	*/
	y = z = -dA;
	x = -dL;
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
			x = -dL;
			z += dA;
		}
		y += dA;
		z = -dA;
	}

	/*
	�������� ��� ������ ��� ���� ������ � ������� ����� �����������
	����� � �������.
	��� �������, ������� �� ������� � ����� ����������, ��������� ���
	��������� ��������.
	*/
	for (int i = 0; i < NP; i++) {

		for (int j = 0; j < 4; j++) {
			particles[i].x[j] = -10.0;
			particles[i].y[j] = -10.0;
			particles[i].z[j] = -10.0;
			particles[i].x_box[j] = 0;
			particles[i].y_box[j] = 0;
			particles[i].z_box[j] = 0;
			particles[i].box_i[j] = 0;
		}

		particles[i].ti = 0;
		particles[i].dt = 0.0;
	}

	/*
	��������� ������ � ���� �������� � "�������" �� ����� ����������
	*/
	for (int i = 0; i < NP; i++) {
		load_information_about_one_particle(loading_file);
	}

	fclose(loading_file);

	/*
	�������������� ������� ������� ������� ��������� � �������������
	��������� ����� ������ �� ��� ������.
	*/
	last = 1;
	time_queue[0].t = 0.0;
	time_queue[0].im = -1;
	for (int i = 1; i < NP; i++) time_queue[i].t = 1.0E+20;

	/*
	������������ ����� ������� ��� ���� ������ � ������������ � ������� "�������"
	*/
	for (int i = 0; i < NP; i++) {
		retime(i);
	}
}


/*
������� ���������� ��������� ������� � ��������� ����.
����� ����������� ������������ ������������� ������ �� �������.

���������:
file_name - ��� �����, � ������� ����� �������� ����������
*/
void save(std::string file_name) {
	double x, y, z, dt;
	double t_global = get_maximum_particle_time();

	FILE *save_file = fopen(file_name.c_str(), "w+");

	// RU: ��������� ���������� � ���������� ������ � �������� �������
	// EN: save information about count of particles and size of the system
	fprintf(save_file, "%d\n", NP);
	fprintf(save_file, "%.15le\n", A);
	fprintf(save_file, "%.15le\n", L);

	// RU: ��������� ���������� ���� ������ � �� ��������: x, y, z, vx, vy, vz
	// EN: we need to save all coordinates of particles: x, y, z, vx, vy, vz
	for (int i = 0; i < NP; ++i) {
		fprintf(save_file, "%d\n", i);

		/*
		�������������� ������� i � ���������� �������� �������,
		� ����� ������ ���� ������ ����� �������� ��� ������ ������� �������,
		������� ��������� �� �������� ���������� ������������� � ������� �������.
		*/
		dt = t_global - particles[i].t;

		/*
		��� ���������� ��������� ������� �� ���������
		���������� ���������� ������ � ��������� �������, ���, �����
		����� ������� ��� � ����� (L/2, A/2), � �� � (0, 0).
		*/
		x = particles[i].x[0] + particles[i].vx * dt;
		y = particles[i].y[0] + particles[i].vy * dt;
		z = particles[i].z[0] + particles[i].vz * dt;
		fprintf(save_file, "%.15le %.15le %.15le\n", x, y, z);
		fprintf(save_file, "%.15le %.15le %.15le\n", particles[i].vx,
			    particles[i].vy, particles[i].vz);
	}

	fclose(save_file);
}


/*
��������� ������ "������" �������.

���������:
NN - ����� ������ � ����� ������ ��������������� ���������, �� ������
�������� �������� ����������� �����
etta - ��������� ���������, ������� ���������� ������ � �������
*/
void new_seed(int NN, double etta) {
	int KZ = 2;

	double axy, axz, v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z;

	// ���������� ����� ����� ������
	double rb = 2.0 * sqrt(2.0);
	double rbp = sqrt(2.0);

	// ������� ���������� ���������� ������
	double A0 = rb * (NN - 1.0) + 2.0;
	double L0 = rb * (2.0 * NN - 1.0) + 2.0;

	double mL0 = L0;

	double Betta = 0.0;

	double XC, YC, ZC; // ���������� ������ ������
	int NA; // ���������� ��������� ������
	double ax, ay, az, vx, vy, vz;

	NP = 2 * NN * (NN * NN + (NN - 1) * (NN - 1)); //  ��������� ����� ������
	NP = NP + (2 * NN - 1) * 2 * (NN - 1) * NN;  //

	// ������� �������������� ���������
	double etta0 = (4.0 * PI * NP) / (3.0 * A0 * A0 * (L0 - 2.0));

	// �������� ������� ������������ ���������� �������
	double Alpha = etta0 / etta;
	double sk = L0 / (2.0 * Alpha * (L0 - 2.0));

	/*
	Betta - ����������� ���������� ������� �� ��������� ��������������������
	�������� � ��������, ������������� �� ����� ������ �������
	*/
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

	/*
	"������������" ��������� ��������� �����, ����� ��� ������ ������� ����������
	���������� ��������������� ����� ���� ������
	*/
	//srand(time(NULL));
	srand(10);

	for (int i = 0; i < 2 * NN; i++)
		for (int j = 0; j < NN; j++)
			for (int k = 0; k < NN / 2; k++) {
				/*
				������ �������� ��������� �������, ��� ����� ���������
				��� ���������� ������ �����(� ��������� �����������),
				����� ��� ����� ���� ������� � ������ �������� �������.
				*/
				vx = double(double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100))
					- double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)));
				vy = double(double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100))
					- double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)));
				vz = double(double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100))
					- double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)));

				for (int ii = -1; ii < 2; ii += 2) {
					if ((i == 0) && (ii > 0)) continue;

					ax = XC + i * ii*rbp + 1.0e-7;

					for (int jj = -1; jj < 2; jj += 2) {
						if ((j == 0) && (jj > 0)) continue;

						ay = YC + j * jj*rbp + 1.0e-7;

						for (int kk = -1; kk < 2; kk += 2) {
							if (((i % 2 == 0) && (j % 2 == 0)) ||
								((i % 2 != 0) && (j % 2 != 0)))
								az = ZC + kk * (rbp + k * rb) + 1.0e-7;
							else {
								if ((k == 0) && (kk > 0)) continue;
								az = ZC + kk * k*rb + 1.0e-7;
							}

							// ��������� ���������� ��� ����� ��������
							particles[NA].x[0] = ax;
							particles[NA].y[0] = ay;
							particles[NA].z[0] = az;

							if ((j == 0) && (fabs(ZC - az) < 0.1E-15L)
								&& ((i / 2) * 2 != i))
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
			particles[ii + jj * NP].x[0] = particles[ii].x[0] + jj*mL0;
			particles[ii + jj * NP].y[0] = particles[ii].y[0];
			particles[ii + jj * NP].z[0] = particles[ii].z[0];

			particles[ii + jj * NP].vx = particles[ii].vx;
			particles[ii + jj * NP].vy = particles[ii].vy;
			particles[ii + jj * NP].vz = particles[ii].vz;
		}
	}

	NP = NP*KZ;

	/*
	"���������" ������� �� ����������� ���������.
	��������� Betta - ��� ����������� ����������,
	�� ���� ���������� ��� ���������� � ��������� ������.
	*/
	for (int ii = 0; ii < NP; ii++) {
		particles[ii].x[0] = particles[ii].x[0]*Betta;
		particles[ii].y[0] = particles[ii].y[0]*Betta;
		particles[ii].z[0] = particles[ii].z[0]*Betta;
	}

	// ������������ � ������� ������ �������� Lx ������� ����� ������:
	double Lx = 0.0;
	for (int i = 0; i < NP; i++) {
		Lx += particles[i].y[0]*particles[i].vz - particles[i].z[0]*particles[i].vy;
	}
	printf("\n Lx = %.15le\n", Lx);

	for (int i = 0; i < NP; i++) {
		particles[i].t = 0.0;
		particles[i].dt = 0.0;
		particles[i].ti = 0;
	}

	// �������� ��� ������� ������� �������� ��������� ���������
	L = ((4.0 * PI * NP) / etta) / (3.0 * A * A) + 2.0;

	/*
	C�������� � ��������� ������� �� ����� ����� ����������������
	��� ����������� ����������.
	*/
	save("new.txt");
	load_seed("new.txt");
}


void clear_cell_from_virt_particle(int i, int j) {
	particle p1 = particles[i];
	Box p1_box = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]];

	if (p1_box.particles[p1.box_i[j]] % NP == i) {
		short end = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].end;

		if (p1.box_i[j] <= end) {
			/*
			������� ������� �� ������ � ������� ��� ����������,
			������ �� ��������� ������� �� ����� ������ ������ � ������ �
			��������� ����� ������ � ������ �� 1.
			*/
			int fj = p1_box.particles[end] % NP;
			int kj = p1_box.particles[end] / NP;

			p1_box.particles[p1.box_i[j]] = p1_box.particles[end];
			particles[fj].box_i[kj] = p1.box_i[j];
			--p1_box.end;
			boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]] = p1_box;
		}
	}
}


/*
������� ��������� ��������� ������ � ������������ � �����������
� ������� ��������.

���������:
im - ����� �������
jm - ����� ������ ������� ��� ����� �������
*/
void reform(int &im, int &jm, int &vp1, int &vp2) {
	particle p1 = particles[im];
	double q1, q2, z, dx, dy, dz;

	/*
	   ���� jm >= 0 �� ����������� ������� - ��� ���������� ����
	   ������ ��� �������.
	*/
	if (jm >= 0) {
		particle p2 = particles[jm];

		for (int j = 0; j < 4; ++j) {
			if (p2.x[j] > 0.0) {
				p2.x[j] += p2.vx * p2.dt;
				p2.y[j] += p2.vy * p2.dt;
				p2.z[j] += p2.vz * p2.dt;
			}
		}
		p2.t = p1.t;
		p2.dt = 0.0;
		dx = p1.x[vp1] - p2.x[vp2];
		dy = p1.y[vp1] - p2.y[vp2];
		dz = p1.z[vp1] - p2.z[vp2];

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

		particles[im] = p1;
		particles[jm] = p2;
	}
	else
	{
		// ���������� � ��������� �������
		if (jm == -1) {
			particles[im].vx = -particles[im].vx;
		}
		else {
			// ������� �������� ����� �������� �������
			if (jm > -10) {
				for (int j = 0; j < 4; ++j) {
					if (p1.x[j] > 0.0) {
						Box p1_box = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]];
						short end = p1_box.end;

						/*
						������� ������� �� ������ � ������� ��� ����������,
						������ �� ��������� ������� �� ����� ������ ������ � ������ �
						��������� ����� ������ � ������ �� 1.
						*/
						int fj = p1_box.particles[end] % NP;
						int kj = p1_box.particles[end] / NP;
						
						/*
						if (p1_box.particles[p1.box_i[j]] % NP != im) {
							printf("\n p1_box.particles[p1.box_i[j]] = %d   p1_box.particles[p1.box_i[j]] / NP %d \n", p1_box.particles[p1.box_i[j]], p1_box.particles[p1.box_i[j]] % NP);
							printf("\n %d %d %d  im %d  j %d\n", p1.y_box[j], p1.z_box[j], p1.x_box[j], im, j);
							throw "EEE";
						}
						*/

						p1_box.particles[p1.box_i[j]] = p1_box.particles[end];
						particles[fj].box_i[kj] = p1.box_i[j];
						--p1_box.end;
						boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]] = p1_box;

						if (jm == -2) {
							--p1.x_box[j];
							if (p1.x[j] > 0.0) {
								p1.x[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].x2;
							}
						}
						if (jm == -4) {
							++p1.x_box[j];
							if (p1.x[j] > 0.0) {
								p1.x[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].x1;
							}
						}

						if (jm == -5) {
							--p1.y_box[j];
							p1.y[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].y2;
						}
						if (jm == -6) {
							++p1.y_box[j];
							p1.y[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].y1;
						}
						if (jm == -7) {
							--p1.z_box[j];
							p1.z[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].z2;
						}
						if (jm == -8) {
							++p1.z_box[j];
							p1.z[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].z1;
						}

						/*
						����� ����, ��� ������ ������ ��� ������� ��� ������,
						��������� ������� � ����� ������ (��� ���������� � � ������
						������, ���� ������ ������ �� ���������).
						*/
						p1.box_i[j] = ++boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].end;
						boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].particles[p1.box_i[j]] = im + NP*j;
						particles[im] = p1;
					}
				}
				if (particles[im].y_box[0] == 0 || particles[im].z_box[0] == 0 ||
					particles[im].y_box[0] == K || particles[im].z_box[0] == K) {
					change_with_virt_particles(im, jm);
				}
			}
			// ����������� ������� Y = A - 2
			if (jm == -11) {
				if (p1.vy > 0.0) {
					create_virt_particle(im, jm);
				}
				else {
					particles[im].x[1] = -10.0;
					particles[im].x[3] = -10.0;

					clear_cell_from_virt_particle(im, 1);
					clear_cell_from_virt_particle(im, 3);
				}
			}
			// ����������� ������� Y = 2 
			if (jm == -12) {
				if (p1.vy < 0.0) {
					create_virt_particle(im, jm);
				}
				else {
					particles[im].x[1] = -10.0;
					particles[im].x[3] = -10.0;

					clear_cell_from_virt_particle(im, 1);
					clear_cell_from_virt_particle(im, 3);
				}
			}
			// ����������� ������� Z = A - 2
			if (jm == -13) {
				if (p1.vz > 0.0) {
					create_virt_particle(im, jm);
				}
				else {
					particles[im].x[2] = -10.0;
					particles[im].x[3] = -10.0;

					clear_cell_from_virt_particle(im, 2);
					clear_cell_from_virt_particle(im, 3);
				}
			}
			// ����������� ������� Z = 2
			if (jm == -14) {
				if (p1.vz < 0.0) {
					create_virt_particle(im, jm);
				}
				else {
					particles[im].x[2] = -10.0;
					particles[im].x[3] = -10.0;

					clear_cell_from_virt_particle(im, 2);
					clear_cell_from_virt_particle(im, 3);
				}
			}
		}
	}
}


void check_overlap() {
	double time = get_maximum_particle_time();

	for (int i = 0; i < particles_for_check_count; i++) {
		for (int j = i + 1; j < particles_for_check_count; j++) {
			particle p01 = particles[particles_for_check[i]];
			particle p02 = particles[particles_for_check[j]];

			double dt01 = time - p01.t;
			double dt02 = time - p02.t;
			double ddx = p01.x[0] + p01.vx*dt01 - p02.x[0] - p02.vx*dt02;
			double ddy = p01.y[0] + p01.vy*dt01 - p02.y[0] - p02.vy*dt02;
			double ddz = p01.z[0] + p01.vz*dt01 - p02.z[0] - p02.vz*dt02;
			double rr = ddx*ddx + ddy*ddy + ddz*ddz;

			if (rr < 4.0 - 1.0e-14) {
				printf("\n r = %.15le ", rr);
				printf("\n dt01 = %.15le < p01.dt = %.15le ? \n ", dt01, p01.dt);
				printf("\n dt02 = %.15le < p02.dt = %.15le ? \n ", dt02, p02.dt);

				throw "ALARM";
			}
		}
	}
}


/*
������� "���", �������� ���� ���������,
���������� ������� �������� ������� � ������� 1 ����������
�� ������ �������  - (NP/2) ������������ � �������.
*/
void step() {
	particle p1;
	int i, im, jm, vp1, vp2;
	double time = get_maximum_particle_time();

	COLL_COUNT = 0;
	jm = 0;

	while (COLL_COUNT < NP / 2 || jm < 0) {
		// ��������� ������ ������� �� ������� �������
		im = time_queue[1].im;
		jm = time_queue[1].jm;
		vp1 = time_queue[1].vp1;
		vp2 = time_queue[1].vp2;

		/*
		int f = 376;
		if (im == f || jm == f) {
			printf("%d %d %d %d \n", im, jm, vp1, vp2);
			printf("particles[%d][0] %.5le %.5le %.5le %d %d %d \n", f, particles[f].x[0], particles[f].y[0], particles[f].z[0], particles[f].x_box[0], particles[f].y_box[0], particles[f].z_box[0]);
			printf("particles[%d][1] %.5le %d %d %d \n", f, particles[f].x[1], particles[f].x_box[1], particles[f].y_box[1], particles[f].z_box[1]);
			printf("particles[%d][2] %.5le %d %d %d \n", f, particles[f].x[2], particles[f].x_box[2], particles[f].y_box[2], particles[f].z_box[2]);
			printf("particles[%d][3] %.5le %d %d %d \n", f, particles[f].x[3], particles[f].x_box[3], particles[f].y_box[3], particles[f].z_box[3]);
			check_particles();
		}
		*/
		
		//printf("%d %d %d %d \n", im, jm, vp1, vp2);
		//check_particles();

		// ������� ������ �������
		delete_event(1);

		p1 = particles[im];

		/*
		�������� ��������� �� ������� ������, ������������ � ����
		�������, ����������� �� -1, ���, ����� ��� �� ��������� ��
		������������ �������.
		*/
		p1.ti = -1;
		if (jm >= 0) particles[jm].ti = -1;

		/*
		���� ��������� ������� -100, ������ ������
		�������������� ������� � ������� ���������
		������ ������� � ���������� �������� �������.
		*/
		if (jm == -100) {
			p1.dt = time - p1.t;
		}

		for (int j = 0; j < 4; ++j) {
			if (p1.x[j] > 0.0) {
				p1.x[j] += p1.vx * p1.dt;
				p1.y[j] += p1.vy * p1.dt;
				p1.z[j] += p1.vz * p1.dt;
			}
		}
		p1.t += p1.dt;
		p1.dt = 0.0;
		time = p1.t;
		
		particles[im] = p1;

		// ���������� ��������� � ������� �������� ������������� �������
		reform(im, jm, vp1, vp2);

		retime(im);

		/*
		���� jm > 0, �� ��� ������, ��� ���������
		���������� ������ � �������, �����������
		������� ����������
		*/
		if (jm >= 0) {
			++COLL_COUNT;
			retime(jm);
		}
	}

	//printf("\n Sync! \n");
	//check_particles();
	/*
	��������� ���������� ������������� ������ �� ������������ �������
	*/
	time = get_minimum_particle_time();
	particle *p = particles;
	for (i = 0; i < NP; ++i, ++p)
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
� ���������� �������� �������.

���������:
steps - ����� ������� ��������� ��� ���������� ������ � ������� ���������
accuracy - ����� ����� �� � �� ������� ����������� ������� ���������
� � ������� ��������� ���������� ������
file_name - ��� ����� ��� ���������� ������
*/
void image(long long int steps, short accuracy, std::string file_name) {
	// ����� �����, � ������� ����� ������������ ���������� ������
	const int W = int(L*accuracy);

	// ������, � ������� ����� �������� ������ �� ���������� ������
	int img[10000];
	double t_global, dt;

	printf("INFO: Image started for %d steps with accuracy %d\n", steps, accuracy);

	// ��������� ������ ������, ������ ��� ������� ���������� ������
	for (short g = 0; g < 10000; g++) img[g] = 0;

	for (short h = 0; h < steps; h++) {
		step();  // ����� ���� ���������� �� ������� ����� ��������

		// ������� ������ � ����� ���������� � ���������
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("%d / %d", int(h + 1), steps);

		/*
		������������ ���������� ����� �������, ������� �����
		������ �������� ������������ ������� ������ � �������.
		*/
		t_global = get_maximum_particle_time();

		for (int i = 0; i < NP; ++i) {
			/*
			�������������� ������� � ���������� �������� �������
			*/
			dt = t_global - particles[i].t;

			/*
			���������� � ����� �� ����� �� X �������� ������ �������
			����� ������������� � �� ������� � ���������� �������� �������
			*/
			int m = int((particles[i].x[0] + particles[i].vx * dt) * 10.0 * accuracy) / 10;
			img[m]++;  // ����������� ���������� ������ � ������ �� 1.
		}
	}

	double x = 0.0, delta_x = 1.0 / accuracy;
	FILE *profile_file = fopen(file_name.c_str(), "w+");
	for (int f = 0; f < W; ++f) {
		/*
		��������� ��������� ���������� X ��� ������ � ����� ������,
		������� ���������� � ������ "����" �� ����� ���������
		*/
		fprintf(profile_file, "%f %d\n", x, img[f]);
		x += delta_x;
	}
	fclose(profile_file);

	printf("\n INFO: Image completed. Information saved to file: %s \n", file_name.c_str());

	save("new.txt");
	load_seed("new.txt");
}


/*
void function_g(int steps_count, double x1, double x2, std::string file_name) {

	int Ni = 0;
	long int g_data[1000];

	for (int i = 0; i < 1000; i++)
		g_data[i] = 0;

	for (int w = 0; w < steps_count; w++) {
		step();

		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("%d / %d", w + 1, steps_count);

		double t_global = get_maximum_particle_time();

		for (int i = 0; i < NP; i++) {
			particle p1 = particles[i];

			double dt = t_global - p1.t;
			p1.x = p1.x + p1.vx * dt;
			p1.y = p1.y + p1.vy * dt;
			p1.z = p1.z + p1.vz * dt;

			if ((L + p1.x >= x1) && (L + p1.x <= x2) &&
				(p1.y < A - 5.0) && (p1.y > 5.0 - A) &&
				(p1.z < A - 5.0) && (p1.z > 5.0 - A))
			{
				Ni += 1;

				for (int j = i + 1; j < NP; j++) {
					particle p2 = particles[j];

					double dt = t_global - p2.t;
					p2.x = p2.x + p2.vx * dt;
					p2.y = p2.y + p2.vy * dt;
					p2.z = p2.z + p2.vz * dt;

					double dx = p1.x - p2.x;
					double dy = p1.y - p2.y;
					double dz = p1.z - p2.z;
					double r = sqrt(dx*dx + dy*dy + dz*dz);

					if (r < 10.0L) {
						int u = int(r * 100);

						g_data[u] = g_data[u] + 1;
					}
				}
			}
		}
	}

	FILE *data_file = fopen(file_name.c_str(), "w+");
	fprintf(data_file, "%d\n", steps_count);

	for (int i = 0; i < 1000; i++)
		fprintf(data_file, "%.5le\n", double(g_data[i]) / ((double(Ni) / steps_count) * 4 * PI * double(i * i / 10000.0)));
	fclose(data_file);

}
*/

/*
void count_of_nearest_particles(int steps_count, double x1, double x2, std::string file_name) {

	double Ni = 0, count_of_particles = 0;

	for (int w = 0; w < steps_count; w++) {
		step();

		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("%d / %d", w + 1, steps_count);

		double t_global = get_maximum_particle_time();

		for (int i = 0; i < NP; i++) {
			particle p1 = particles[i];

			double dt = t_global - p1.t;
			p1.x = p1.x + p1.vx * dt;
			p1.y = p1.y + p1.vy * dt;
			p1.z = p1.z + p1.vz * dt;

			if ((L + p1.x >= x1) && (L + p1.x <= x2) &&
				(p1.y < A - 5.0) && (p1.y > 5.0 - A) &&
				(p1.z < A - 5.0) && (p1.z > 5.0 - A))
			{
				Ni += 1;

				for (int j = 0; j < NP; j++) {
					if (i == j) continue;

					particle p2 = particles[j];

					double dt = t_global - p2.t;
					p2.x = p2.x + p2.vx * dt;
					p2.y = p2.y + p2.vy * dt;
					p2.z = p2.z + p2.vz * dt;

					double dx = p1.x - p2.x;
					double dy = p1.y - p2.y;
					double dz = p1.z - p2.z;
					double r = sqrt(dx*dx + dy*dy + dz*dz);

					if (r < 3.0L) {
						count_of_particles += 1;
					}
				}
			}
		}
	}

	FILE *data_file = fopen(file_name.c_str(), "w+");
	fprintf(data_file, "%.15le\n", double(count_of_particles / Ni));
	fclose(data_file);

}
*/

/*
������� �������� ������� ������� � ����,
��������� ������������� ������������ ������ � ����,
������������ ��������� ������.

���������:
x1 - ��������� X ���������� "�����"
x2 - �������� X ���������� "�����"
dots_per_particle - ���������� ������ ������ ����� ������ ���������� ���������� � �������
steps - ����� ���������� ����� ����� �������� ������ � ��������� ������ � ��������� ����
file_name - ��� ����� ��� ���������� ������
*/
void profile(double x1, double x2, int dots_per_particle, long long int steps, std::string file_name) {
	double x, y, z, dt;
	int i, j;

	// ��������� ���� �� ������
	FILE *profile_file = fopen(file_name.c_str(), "w+");

	printf("INFO: Profile started.\n");

	for (i = 0; i < dots_per_particle; ++i) {
		/*
		����� ����������� ������ � ���������� ������
		������ ��������� ���������� ���������� �� ������� � �������
		*/
		for (j = 0; j < steps; j++)
			step();

		double t_global = get_maximum_particle_time();

		for (j = 0; j < NP; j++) {
			/*
			�������������� ������� �� ������� � ����������
			�������� ������� � ��������� � ���������� �,
			���� ������� ��������� � ��������� ���������,
			�� ��������� � Y � Z ���������� � ����.
			*/
			dt = t_global - particles[j].t;
			x = particles[j].x[0] + particles[j].vx * dt;

			if (x <= x2 && x >= x1) {
				y = particles[j].y[0] + particles[j].vy * dt;
				z = particles[j].z[0] + particles[j].vz * dt;

				// ��������� Y � Z ���������� �������
				fprintf(profile_file, "%.15le\n", y);
				fprintf(profile_file, "%.15le\n", z);
			}
		}
	}

	fclose(profile_file);

	printf("INFO: Profile completed. Information saved to file: %s\n", file_name.c_str());

	save("new.txt");
	load_seed("new.txt");
}


/*
������� ������ �������, ��������� ������� ��� ��������� ������� �� �������� ���������
� �������� ������������ ����� �� ���������.

���������:
compress_to_etta - �������� ���������
delta_L - ������������ ����������� ��� �� L
steps - ���������� ���������� �� ���� ������� � ������� �����
���������� �������� �� ���������
type - ��� ������:
0 - ����������� ������ ����� ������
1 - ����������� ������ ������ ������
2 - ����������� ������������ ��� ������ �� ���������� ����������
*/

void compress(double compress_to_etta, double delta_L, int KK1, int KK2, int KK3, int type) {
	int compression_steps_done = 0;
	double etta = (4.0 * PI * NP) / (3.0 * A * A * (L - 2.0));
	double min = 1.0e+100, max = -1.0, x, dL, dx;
	double t_global, dt;
	double L_ideal = ((4.0 * PI * NP) / compress_to_etta) / (3.0 * A * A) + 2.0;

	printf("\n INFO: Start to change system density... \n");

	printf(" \n  Maximum delta_L = %.15le \n", delta_L);
	printf("  Program will wait %d collissions per particle between each changes in density. \n", KK1);
	printf("  Program will wait %d collissions per particle between each %d changes in density. \n", KK3, KK2);
	printf("  Type of compression: ");
	if (type == 0) printf("only position of left wall will be changed");
	if (type == 1) printf("only position of right wall will be changed");
	if (type == 2) printf("position of both walls will be changed \n\n");

	printf("etta = %.15le, should be equal to %.15le \n", etta, compress_to_etta);

	// ����� ��������� � ��������� � 12 ������
	while (fabs(etta - compress_to_etta) > 1.0e-12) {
		if (etta < compress_to_etta) {
			max = -100.0;
			min = 1.1e+10;

			
			//������� ��� ������� � ������� ����� ������� ����� ��������������
			//��� ������� � ����� ������� ���������� ������, �������� ������
			//������������� �� ��������� ������, ����� ����� �� ������� �����
			//����������� ������ �� ���������� ������.
			
			t_global = get_maximum_particle_time();

			for (int i = 0; i < NP; ++i) {
				dt = t_global - particles[i].t;
				x = particles[i].x[0] + particles[i].vx * dt;
				if (x < min) min = x;
				if (x > max) max = x;
			}

			
			//������������ �� ������� ������ ������� ��������� � ������
			//� ������ ��������� �������
			
			min = min - 1.0;
			max = L - max - 1.0;

			// �������� ���������� � ���� ���������� � ��������� � dL
			dL = min;
			if (dL > max) dL = max;

			// ������� �� ������� � �������� � �� ������� ������
			dL = dL / 1.1;
			if (dL < 0.1e-12) dL = 0.01e-30; // �� ������� ������ ���� ������� ������
			if (dL > delta_L) dL = delta_L;  // �������� ������ �� ������ ��� �� delta_L

			dx = dL / 2.0;  // ������������ �������� ��� ���� ������

			
			//���� �������� �������� ������ ����� ������� ������
			//��� ��������� (�.�. ������������� ��������� ������ � ������� etta
			//����� ������ / ������ ���, ������� ���� ������� � ����������),
			//�� ���������� ������ ������ �������� L, ����� ��������� �������
			//� ��������� ����������.

			if (L - dL < L_ideal) {
				dL = 0.0;
				L = L_ideal;
				dx = 0.0;
			}
		}
		else {
			
			//���� ������� ���������� ���������, �� ������ ����
			//���������� ����������� ��� ��� ��������� ��������� ������
			//� �������� ��������� ������ ��� ��� ������.
			
			dL = L - L_ideal;
			if (dL < -delta_L) dL = -delta_L; // ��� ���������� �������
			dx = dL / 2.0;
		}

		if (type != 1) {
			
			//������� ��� ������� �����, ����� ������� ���������� ������ ����� ������.
			//���������� ������ �� ������ ������ �� ���������, �.�. ����� ����������� ������
			//����� �� �������� ��� ������ �� �� �� ���������� - � ����� �� ������� ���
			//������� ����� �� dL/2.0.
			
			for (int w = 0; w < NP; ++w) {
				for (int j = 0; j < 4; j++) {
					particles[w].x[j] -= dx;
				}
			}

			if (type == 0) {
				dL = dx;
			}
		}

		// �������� ������� - ���������� ���������� ��������� ��������� ���� ������
		L -= dL;

		
		//��������� ��������� ������� � ����� ��������� ��� ���������� �����
		//��������� � ������� ����������� �������������
		
		save("tmp");
		load_seed("tmp");

		
		//����� ������� �������� ������ ������ ��������� ���������� ����������
		//�� ������� � �������
		for (short i = 0; i < KK1; i++)
			step();

		compression_steps_done += 1;
		if (compression_steps_done % KK2 == 0) {
			for (short i = 0; i < KK3; i++)
				step();
		}

		// ������������ ��������� ����� ������
		etta = (4.0 * PI * NP) / (3.0 * A * A * (L - 2.0));

		printf("etta = %.15le, should be equal to %.15le \n", etta, compress_to_etta);
	}

	printf("\n INFO: System density was sucessfully changed to %.15le \n", etta);
}

/*
������� ���������� ������������ �� ���������� �
��������� ����� ��������.

���������:
file_name - ��� ����� � ��������� ����� ������������ ��� ���������.
*/
void init(std::string file_name) {
	using namespace std;
	clock_t start, end, result;
	char command[255], parameter[255];
	long long int i, steps;
	ifstream command_file(file_name.c_str());

	while (!command_file.eof()) {
		command_file.getline(command, 255, ' ');
		string str_command = command;

		printf("\n\n<==========================>\n");

		// ���� ���������� ������� ����� �����
		if (str_command.compare("new") == 0) {
			int NN;
			double etta;

			command_file >> NN;
			command_file >> etta;
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������
			new_seed(NN, etta);

			print_system_parameters();
		}
		// ���� ���������� ��������� ��������� ������� �� �����
		if (str_command.compare("load") == 0) {
			command_file.getline(parameter, 255, '\n');
			load_seed(parameter);

			printf("\n System was successfully loaded from file '%s' \n", parameter);

			print_system_parameters();
		}
		// ���� ���������� ���������� �������� ������� � �������
		// ���������� ���������� ����������
		if (str_command.compare("step") == 0) {
			command_file >> steps;
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������

			printf(" INFO: Step start. \n");

			start = clock();

			for (i = 0; i < steps; ++i) {
				step();
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
				printf("%d / %d", int(i + 1), steps);

				// ���� �� ����� ���� � �������� �������,
				// �� ������� ��� ������ 1000 ���������� �� �������,
				// ����� �� �� �������� ���������� ��������.
				if (i % 1000 == 0) {
					FILE *history_file = fopen("history.txt", "w+");
					fclose(history_file);

					check_particles();
				}
			}

			end = clock();
			result = end - start;

			printf("\n INFO: finished %d collisions per particle \n", steps);
			printf("Total Time = %f seconds. \n", double(result / CLOCKS_PER_SEC));

			/*
			������������ ������ �������� ������� Lx � ������� ��� �� �����.
			*/
			double Lx = 0.0;
			for (int i = 0; i < NP; i++) {
				Lx += particles[i].y[0]*particles[i].vz - particles[i].z[0]*particles[i].vy;
			}
			printf("\n Lx = %.15le\n", Lx);
		}
		// ���� ���������� ������� ������ �� ������� ��������� � �������
		if (str_command.compare("image") == 0) {
			command_file >> steps; // ����� ���������� �� ���� ������� �� ����� ���������
			command_file >> i; // ��������. ����� ����� ������� �� ���� ������ �������
			command_file.getline(parameter, 255, '\n');
			image(steps, i, parameter);
		}
		// ���� ���������� ������� ������ �� "�������" � �������
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
		// ���� ���������� ��������� ��������� �������
		if (str_command.compare("save") == 0) {
			command_file.getline(parameter, 255, '\n');
			save(parameter);
			printf("\n INFO: particles coordinates saved to '%s' \n", parameter);
		}
		if (str_command.compare("function_g") == 0) {
			double x1, x2;
			command_file >> steps;
			command_file >> x1;
			command_file >> x2;
			command_file.getline(parameter, 255, '\n');
			//function_g(steps, x1, x2, parameter);
			printf("\n Function g finished \n");
		}
		if (str_command.compare("count_of_nearest_particles") == 0) {
			double x1, x2;
			command_file >> steps;
			command_file >> x1;
			command_file >> x2;
			command_file.getline(parameter, 255, '\n');
			//count_of_nearest_particles(steps, x1, x2, parameter);
			printf("\n Function count_of_nearest_particles finished \n");
		}
		// ���� ���������� �������� ��������� �������
		if (str_command.compare("compress") == 0) {
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������
			printf("\n Sorry, 'compress' command was deprecated, we need to use \n");
			printf("'compress_two_walls'(short form: 'compresst') instead. You can also use ");
			printf("'compress_left_wall'(short form: 'compressl') \n or 'compress_right_wall'");
			printf("(short form: 'compressr') commands to change system density.\n");
			exit(1);
		}
		// ���� ���������� �������� ��������� �������
		// �������, ������� ��� ��������� ������
		if ((str_command.compare("compress_two_walls") == 0) ||
			(str_command.compare("compresst") == 0)) {
			int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // ��������� ���������
			command_file >> delta_L;  // ����������� ���������� �������� ����������� ������
			command_file >> KK1;  // ���������� ���������� ����� ������� ���� ������
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������
			compress(etta, delta_L, KK1, KK2, KK3, 2);

			print_system_parameters();
		}
		// ���� ���������� �������� ��������� �������
		// ������� ������ "�����" ��������� �����, x = 0
		if ((str_command.compare("compress_left_wall") == 0) ||
			(str_command.compare("compressl") == 0)) {
			int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // ��������� ���������
			command_file >> delta_L;  // ����������� ���������� �������� ��������� ���������
			command_file >> KK1;  // ���������� ���������� ����� ������� ���� ������
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������
			compress(etta, delta_L, KK1, KK2, KK3, 0);

			print_system_parameters();
		}
		// ���� ���������� �������� ��������� �������
		// ������� ������ "������" ��������� �����, x = L
		if ((str_command.compare("compress_right_wall") == 0) ||
			(str_command.compare("compressr") == 0)) {
			int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // ��������� ���������
			command_file >> delta_L;  // ����������� ���������� �������� ��������� ���������
			command_file >> KK1;  // ���������� ���������� ����� ������� ���� ������
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(parameter, 255, '\n');  // ��������� ���������� ������
			compress(etta, delta_L, KK1, KK2, KK3, 1);

			print_system_parameters();
		}
		if (str_command.empty()) break;
	}
	FILE *result_flag = fopen("result", "w+");
	fclose(result_flag);
}

/*
��������� ������� ��� ���������, ������� ���� "program.txt"
� ������ ��������� �����������, ��������� � ���� �����.
*/
int main()
{
	FILE *history_file = fopen("history.txt", "w+");
	fclose(history_file);

	particles_for_check[0] = 3026;
	particles_for_check[1] = 6974;
	particles_for_check[2] = 13500;
	particles_for_check_count = 0;

	init("program.txt");
	return 0;
}