/*

Copyright 2006-2017 Vladimir Veshnev, Aleksey Geraskin, Marina Nurlygayanova,
Timur Nurlygayanov

This file is part of MD program.

MD is free software : you can redistribute it and / or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MD. If not, see <http://www.gnu.org/licenses/>

*/

/*
ÈÌÏÎÐÒÈÐÓÅÌÛÅ ÌÎÄÓËÈ
*/
#include "stdafx.h"  // ñòàíäàðòíàÿ áèáëèîòåêà
#include <time.h>    // áèáëèîòåêà äëÿ ïîäñò÷¸òà âðåìåíè
#include <fstream>   // áèáëèîòåêà äëÿ ðàáîòû ñ ôàéëàìè
#include <cmath>     // áèáëèîòåêà äëÿ âûçîâà ìàòåìàòè÷åñêèõ ôóíêöèé


/*
ÎÁÚßÂËÅÍÈÅ ÏÅÐÅÌÅÍÛÕ
*/

// Êîëè÷åñòâî ÿ÷ååê ïî Y, Z è X(K2), íà êîòîðûå ðàçáèâàåòñÿ îáú¸ì.
// Êîëè÷åñòâî ÿ÷ååê ðàñ÷èòûâàåòñÿ äèíàìè÷åñêè ïðè çàãðóçêå ñèñòåìû
// èëè ïðè ïîñåâå (ñì. ôóíêöèþ load).
long int K, K2;

// ÷èñëî ÷àñòèö âî âñ¸ì îáú¸ìå, ïåðåîïðåäåëÿåòñÿ â ôóíêöèÿõ load è new_seed
long int NP;

// ãëîáàëüíûé ñ÷¸ò÷èê ñòîëêíîâåíèé â ñèñòåìå
long int COLL_COUNT = 0;

bool DIFFUSE_RIGHT_WALL = false;

#define PI 3.141592653589793238462
double particle_R = 0.5;
double particle_Rx2 = 2.0*particle_R;
double particle_D2 = 4.0*particle_R*particle_R;
double particle_R3 = particle_R*particle_R*particle_R;
double x_diffuse = 3.0;

bool temp_save;
char temp_save_file[255];

// ïàðàìåòðû îáú¸ìà, çàäàþòñÿ â load()
double A, dA, L, dL;

// Ãëîáàëüíàÿ ïåðåìåííàÿ äëÿ ïîäñ÷¸òà îáùåé êèíåòè÷åñêîé ýíåðãèè âñåõ ÷àñòèö
double global_E = 0.0;

double dV_left_wall, dV_right_wall, N_collissions_left_wall, N_collissions_right_wall, dP_time;

// èíäåêñ ïîñëåäíåãî ýëåìåíòà â î÷åðåäè ñîáûòèé.
long int last;

// îáúåêò "ñîáûòèå"
typedef struct Event_ {
	double t;
	long int im, jm, vp1, vp2;
} Event;

// î÷åðåäü ñîáûòèé - îïòèìàëüíî 16384 ýëåìåíòà
// (ýòî äîëæíî áûòü ÷èñëî-ñòåïåíü äâîéêè, áîëüøåå ÷åì ìàêñèìàëüíîå ÷èñëî ÷àñòèö)
//Event *time_queue = new Event[262144];
Event *time_queue = new Event[16384];

// îáúåêò "÷àñòèöà"
// x, y, z - êîîðäèíàòû ÷àñòèöû
// vx, vy, vz - ïðîåêöèè ñêîðîñòåé ÷àñòèöû
// t - ñîáñòâåííîå âðåìÿ ÷àñòèöû
// dt - âðåìÿ äî áëèæàéøåãî ñîáûòèÿ ýòîé ÷àñòèöû
// max_x - ìàêñèìàëüíàÿ êîîðäèíàòà Õ äî êîòîðîé äîëåòàåò
//         ÷àñòèöà ïðè äèôôóçíîì îòðàæåíèè
// x_box, y_box, z_box - íîìåð ÿ÷åéêè, â êîòîðîé íàõîäèòñÿ ÷àñòèöà
// ti - íîìåð ñîáûòèÿ ÷àñòèöû â äåðåâå ñîáûòèé
// box_i - íîìåð ÷àñòèöû â ÿ÷åéêå
typedef struct particle_ {
	double x[4], y[4], z[4], vx, vy, vz, t, dt, x_max;
	long int x_box[4], y_box[4], z_box[4], ti, box_i[4];
} particle;

// ìàññèâ ÷àñòèö, ðàçìåð ìàññèâà N + îêðóãëåíèå
// â áîëüøóþ ñòîðîíó ê ÷èñëó äàþùåå ñòåïåíü äâîéêè.
//particle *particles = new particle[262144];
particle *particles = new particle[16384];

// êëåòêà. Îáú¸ì ñèñòåìû ðàçäåë¸í íà ìíîæåñòâî êëåòîê,
// êàæäàÿ êëåòêà ñîäåðæèò â ñåáå íåñêîëüêî âèðòóàëüíûõ ÷àñòèö
// x1, y1, z1, x2, y2, z2 - êîîðäèíàòû êîíöà è íà÷àëà êàæäîé ÿ÷åéêè
// particles[] - ñïèñîê âñåõ ÷àñòèö, íàõîäÿùèõñÿ â äàííîé ÿ÷åéêå
// end - èíäåêñ ïîñëåäíåé ÷àñòèöû â ñïèñêå ÷àñòèö äàííîé ÿ÷åéêè
typedef struct Box_ {
	double x1, y1, z1, x2, y2, z2;
	long int particles[32];
	int end;
} Box;

// ìàññèâ êëåòîê äëÿ âñåãî îáú¸ìà
Box boxes_yz[100][100][1000];


// ìàññèâ ñ íîìåðàìè ÷àñòèö, äëÿ êîòîðûõ íàäî ñîõðàíÿòü èñòîðèþ ñîáûòèé
// èñïîëüçóåòñÿ íà ñëó÷àé îòëàäêè ïðîãðàììû äëÿ ñîõðàíåíèÿ èñòîðèè ñîáûòèé
// âûáðàííûõ ÷àñòèö
long int particles_for_check[100];
long int particles_for_check_count = 0;

double *angles = new double[200];

/*
Ýòà ôóíêöèÿ âûâîäèò íà ýêðàí ïàðàìåòðû òåêóùåé ñèñòåìû:
A - ðàçìåð ñèñòåìû (îáëàñòè îáú¸ìà, ãäå ìîæåò íàõîäèòüñÿ öåíòð ÷àñòèöû) ïî Y è Z.
L
N - ÷èñëî ÷àñòèö â ñèñòåìå
etta - ñðåäíÿÿ îòíîñèòåëüíàÿ ïëîòíîñòü ñèñòåìû
*/
void print_system_parameters() {
	double etta = (4.0 * PI * particle_R3 * NP) / (3.0 * A * A * (L - 2.0*particle_R));
	printf("\n\n| A    | %.15le \n| L    | %.15le \n", A, L);
	printf("| N    | %ld                \n", NP);
	printf("| R    | %.15le \n", particle_R);
	printf("| etta | %.15le \n", etta);
}


void change_R(double newR) {
	particle_R = newR;
	particle_Rx2 = 2.0*particle_R;
	particle_D2 = 4.0*particle_R*particle_R;
	particle_R3 = particle_R*particle_R*particle_R;
}


/*
Ôóíêöèÿ âîçâðàùàåò íàèáîëüøåå ñîáñòâåííîå âðåìÿ ÷àñòèö â ñèñòåìå,
÷òî ïîçâîëÿåò ñèíõðîíèçîâàòü âñå ÷àñòèöû ïî âðåìåíè
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
Ôóíêöèÿ äëÿ ïðîâåðêè ñîñòîÿíèÿ ñèñòåìû, â íåé ìû ïðîâåðÿåì, ÷òî
âñå ÷àñòèöû íàõîäÿòñÿ âíóòðè ñèñòåìû, äëÿ êàæäîé ÷àñòèöû ó íàñ ðàññ÷èòàíî
áëèæàéøåå ñîáûòèå, ÷àñòèöû íàõîäÿòñÿ â ïðàâèëüíûõ ÿ÷åéêàõ ñèñòåìû è ïð.

Â ñëó÷àå âîçíèêíîâåíèÿ ëþáûõ ïðîáëåì äàííàÿ ôóíêöèÿ ïðåêðàùàåò ðàáîòó ïðîãðàììû
è âûâîäèò äîïîëíèòåëüíóþ èíôîðìàöèþ îá îáíàðóæåííîé ïðîáëåìå.

Ôóíêöèÿ ðàáîòàåò ìåäëåííî, íåîáõîäèìî èñïîëüçîâàòü â öåëÿõ ïðîâåðêè
èçìåíåíèé â ïðîãðàììå.
*/
int check_particles() {
	double E = 0.0;
	double dt = 0.0, dt2 = 0.0;
	double t_global = get_maximum_particle_time();

	for (long i = 0; i < K + 1; i++) {
		for (long j = 0; j < K + 1; j++) {
			for (long f = 0; f < K2 + 1; f++) {
				if (boxes_yz[i][j][f].end > 30) {
					printf("!! %ld %ld %ld  end = %d", i, j, f, boxes_yz[i][j][f].end);
					throw "AAA";
				}

				for (long int p = 0; p < boxes_yz[i][j][f].end + 1; p++) {
					if (boxes_yz[i][j][f].particles[p] > NP * 4) {
						printf("\n !! %ld %ld %ld  boxes_yz[i][j][f].particles[p] %ld \n", i, j, f, boxes_yz[i][j][f].particles[p]);
						throw "AAA";
					}
				}
			}
		}
	}

	for (long int i = 0; i < NP; i++) {
		particle p1 = particles[i];
		Box p1_box = boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]];

		dt = t_global - p1.t;
		p1.x[0] += p1.vx * dt;
		p1.y[0] += p1.vy * dt;
		p1.z[0] += p1.vz * dt;

		E += p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;

		particle p0 = particles[i];
		for (long j = i + 1; j < NP; j++) {
			particle p2 = particles[j];

			dt2 = t_global - p2.t;

			for (long h1 = 0; h1 < 4; h1++) {
				if (p0.x[h1] < 0.0) continue;
				for (long h2 = 0; h2 < 4; h2++) {
					if (p2.x[h2] < 0.0) continue;

					double dx = p0.x[h1] + p0.vx * dt - p2.x[h2] - p2.vx * dt2;
					double dy = p0.y[h1] + p0.vy * dt - p2.y[h2] - p2.vy * dt2;
					double dz = p0.z[h1] + p0.vz * dt - p2.z[h2] - p2.vz * dt2;
					double r = dx*dx + dy*dy + dz*dz;

					if ((r < particle_D2) && (particle_D2 - 1.0e-5 - r > 0.0)) {
						printf("\n\n particles %ld[%ld] %ld[%ld] r = %.15le", i, h1, j, h2, r);
						throw "Particles overlaps!";
					}
				}
			}
		}

		if (p1.x[1] < 0.0 && (p1.y[0] < 2.0*particle_R || p1.y[0] > A - 2.0*particle_R)) {
			printf("\n\n %ld p1.y[0] %.15le", i, p1.y[0]);
			throw "aa";
		}
		if (p1.x[2] < 0.0 && (p1.z[0] < 2.0*particle_R || p1.z[0] > A - 2.0*particle_R)) {
			printf("\n\n %ld p1.z[0] %.15le", i, p1.z[0]);
			throw "aa";
		}

		// ïðîâåðÿåì èíäåêñû ÿ÷ååê äëÿ âñåõ ÷àñòèö
		if ((p1.x_box[0] > K2) || (p1.y_box[0] > K) || (p1.z_box[0] > K) ||
			(p1.x_box[0] < 0) || (p1.y_box[0] < 0) || (p1.z_box[0] < 0)) {
			throw "Particle locates in incorrect cell.";
		}

		// ïðîâåðÿåì ÷òî ÷àñòèöû íàõîäÿòñÿ â îáú¸ìå
		if ((p1.x[0] > L && p1.x[0] < 0.0) ||
			(p1.y[0] > A && p1.y[0] < 0.0) ||
			(p1.z[0] > A && p1.z[0] < 0.0)) {
			printf(" \n Particle %ld, %.15le, %.15le, %.15le \n ", i, p1.x[0], p1.y[0], p1.z[0]);
			throw "Particle is out of the system boundaries.";
		}

		// ïðîâåðÿåì ÷òî ÷àñòèöû íàõîäÿòñÿ â ïðàâèëüíûõ ÿ÷åéêàõ
		for (int j = 0; j < 4; j++) {
			if (particles[i].x[j] > 0.0) {

				Box p_box = boxes_yz[particles[i].y_box[j]][particles[i].z_box[j]][particles[i].x_box[j]];

				if (((p1.x[j] < p_box.x1) && (abs(p1.x[j] - p_box.x1) > 1.0e-14)) ||
					((p1.x[j] > p_box.x2) && (abs(p1.x[j] - p_box.x2) > 1.0e-14)) ||
					((p1.y[j] < p_box.y1) && (abs(p1.y[j] - p_box.y1) > 1.0e-14)) ||
					((p1.y[j] > p_box.y2) && (abs(p1.y[j] - p_box.y2) > 1.0e-14)) ||
					((p1.z[j] < p_box.z1) && (abs(p1.z[j] - p_box.z1) > 1.0e-14)) ||
					((p1.z[j] > p_box.z2) && (abs(p1.z[j] - p_box.z2) > 1.0e-14))) {

					printf("Particle is out of the box %ld[%d] \n", i, j);

					printf("Boundaries:\n");
					printf("X : [%.15le ; %.15le]\n", p_box.x1, p_box.x2);
					printf("Y : [%.15le ; %.15le]\n", p_box.y1, p_box.y2);
					printf("Z : [%.15le ; %.15le]\n", p_box.z1, p_box.z2);
					printf("x, y, z: %.15le %.15le %.15le\n", p1.x[j], p1.y[j], p1.z[j]);

					printf("p1.t = %.15le, p.im = %ld, p1.jm = %ld \n", p1.t, time_queue[p1.ti].im, time_queue[p1.ti].jm);

					throw "Particle is out of the cell boundary.";
				}
			}
		}

		// Ïðîâåðÿåì ÷òî ÷àñòèöà çàïèñàíà â îäíóþ èç ÿ÷ååê â ñèñòåìå
		bool w = false;
		for (long int ty = 0; ty <= boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].end; ++ty) {
			if (boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].particles[ty] == i) w = true;
		}
		if (w == false) {
			for (long int t = 0; t <= boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].end; ++t) {
				printf("\n %ld ", boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].particles[t]);
			}

			printf("\n Particle %ld doesn't store in the cell.", i);
			printf("\n x %.5le, y %.5le, z %.5le", particles[i].x[0], particles[i].y[0], particles[i].z[0]);
			printf("\n x_box %ld, y_box %ld, z_box %ld", particles[i].x_box[0], particles[i].y_box[0], particles[i].z_box[0]);

			throw "Particle doesn't store in the cell.";
		}

		/*
		// Ïðîâåðÿåì íà ïðàâèëüíîå ëè ñîáûòèå â ëèíåéêå âðåì¸í ññûëàåòñÿ ÷àñòèöà
		if (time_queue[p1.ti].im != i && time_queue[p1.ti].jm != i) {
		Event e = time_queue[p1.ti];
		printf("\n i = %d ; im = %d ; jm = %d ; ti = %d ", i, e.im, e.jm, p1.ti);
		throw "Particle has no correct link to the event.";
		}
		*/
	}

	// Ïðîâåðÿåì òåêóùåå çíà÷åíèå ãëîáàëüíîé êèíåòè÷åñêîé ýíåðãèè ñèñòåìû ñî
	// çíà÷åíèåì ýíåðãèè, êîòîðîå áûëî ïðè çàãðóçêå ñèñòåìû èç ôàéëà â load
	if (abs(E - global_E) > global_E*0.1e-8) {
		printf("\nENERGY was changed: \n E_seed = %.15le \n E_now= %.15le \n", global_E, E);
		throw "ENERGY was changed.";
	}

	return 0;
}


/*
Ôóíêöèÿ ïîäú¸ìà ýëåìåíòà ïî î÷åðåäè ñîáûòèé ê íà÷àëó î÷åðåäè

Àðãóìåíòû:
i - ïîçèöèÿ, íà êîòîðîé íàõîäèòñÿ ýëåìåíò â äàííûé ìîìåíò
t - âðåìÿ äî íàñòóïëåíèÿ äàííîãî ñîáûòèÿ
*/
int get_up(long int i, double &t) {
	long int j = i >> 1;  // ýòî ñäâèã âïðàâî, òî æå ñàìîå ÷òî j = i/2, òîëüêî áûñòðåå
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
Ôóíêöèÿ äîáàâëåíèÿ ñîáûòèÿ â î÷åðåäü ñîáûòèé

Àðãóìåíòû:
i - ÷àñòèöà, ñ ó÷àñòèåì êîòîðîé ïðîèçîéä¸ò íîâîå ñîáûòèå
j - íîìåð ÷àñòèöû ñ êîòîðîé ñòîëêí¸òñÿ ÷àñòèöà i èëè íîìåð ñîáûòèÿ
ñîóäàðåíèÿ ñî ñòåíêîé èëè ïðîõîæäåíèÿ ÷åðåç ïåðèîäè÷åñêèå
ãðàíè÷íûå óñëîâèÿ èëè ìåæäó ÿ÷åéêàìè ñèñòåìû
*/
void add_event(long int &i, long int &j, long int &vp1, long int &vp2) {
	/*
	Ðàññ÷èòûâàåì ïîëíîå âðåìÿ íîâîãî ñîáûòèÿ îò íà÷àëà îòñ÷¸òà
	ãëîáàëüíîãî âðåìåíè ñèñòåìû, òàêèì îáðàçîì ìû ïîëó÷àåì âðåìÿ t,
	ñðàâíèâàÿ êîòîðîå ìû ìîæåì îïðåäåëèòü êàêîå èç ñîáûòèé â ñèñòåìå
	ïðîèçîéä¸ò ðàíüøå
	*/
	double t = particles[i].t + particles[i].dt;

	/*
	Íàõîäèì ïîçèöèþ â äåðåâå âðåì¸í äëÿ íîâîãî ñîáûòèÿ.
	Èçíà÷àëüíî ïîìåùàåì ýòî ñîáûòèå âíèç äåðåâà è ïîçâîëÿåì
	åìó ïîäíÿòüñÿ ïî äåðåâó, åñëè äàííîå ñîáûòèå ïðîèçîéä¸ò ðàíüøå ÷åì
	äðóãèå ñîáûòèÿ
	*/
	particles[i].ti = get_up(last, t);

	/*
	Åñëè íîâîå ñîáûòèå - ýòî ñîáûòèå ñòîëêíîâåíèÿ äâóõ ÷àñòèö, òî
	íåîáõîäèìî äëÿ âòîðîé ÷àñòèöû ñîõðàíèòü äàííûå î å¸ íîâîì ñîáûòèè
	*/
	if (j >= 0) {
		particles[j].dt = particles[i].dt;
		particles[j].ti = particles[i].ti;
	}

	// çàïèñûâàåì íîâîå ñîáûòèå â âûáðàííóþ ÿ÷åéêó â äåðåâå âðåì¸í
	time_queue[particles[i].ti].im = i;
	time_queue[particles[i].ti].jm = j;
	time_queue[particles[i].ti].vp1 = vp1;
	time_queue[particles[i].ti].vp2 = vp2;
	time_queue[particles[i].ti].t = t;

	// óâåëè÷èâàåì ÷èñëî ñîáûòèé íà 1
	last++;
}


/*
Ôóíêöèÿ óäàëåíèÿ ñîáûòèÿ èç î÷åðåäè ñîáûòèé

Àðãóìåíòû:
i - ïîçèöèÿ óäàëÿåìîãî ýëåìåíòà â î÷åðåäè ñîáûòèé
*/
void delete_event(long int i) {
	long int j = i << 1;
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
Ôóíêöèÿ óäàëåíèÿ ñîáûòèé ÷àñòèöû èç ëèíåéêè ñîáûòèé

Àðãóìåíòû:
i - íîìåð ÷àñòèöû èëè îáðàçà, ñîáûòèÿ êîòîðîãî íåîáõîäèìî óäàëèòü
*/
void clear_particle_events(long int &i) {
	// particles[i].ti - ññûëêà íà èíäåêñ ñîáûòèÿ, êîòîðîå õðàíèò ñàìà ÷àñòèöà
	int f = particles[i].ti;

	if (f > 0) {
		long int e = -100, pv1 = 0, pv2 = -1;
		long int kim = time_queue[f].im;
		long int kjm = time_queue[f].jm;

		if (time_queue[f].im == i) {
			/*
			Åñëè ìû óäàëÿåì ñîáûòèå ñòîëêíîâåíèÿ äâóõ ÷àñòèö, òî
			íåîáõîäèìî ïåðåìåñòèòü âòîðóþ ÷àñòèöó â òî æå âðåìÿ,
			â êîòîðîì íàõîäèòñÿ ÷àñòèöà, äëÿ êîòîðîé ìû óäàëÿåì
			ýòî ñîáûòèå è çàìåíèòü ñîáûòèå ñòîëêíîâåíèÿ äâóõ ÷àñòèö
			íà ñîáûòèå "-100", êîãäà âòîðàÿ ÷àñòèöà ïðîñòî äîëåòèò äî
			ìåñòà ïðåäïîëàãàåìîãî ñòîëêíîâåíèÿ è ïîñëå ýòîãî äëÿ íå¸
			áóäåò ðàññ÷èòàíî íîâîå ñîáûòèå.
			*/
			if (time_queue[f].jm >= 0) {
				double dt = particles[kim].t - particles[kjm].t;  // ðàçíèöà â òåêóùåì âðåìåíè

				for (long int r = 0; r < 4; ++r) {
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
			double dt = particles[kjm].t - particles[kim].t;  // ðàçíèöà â òåêóùåì âðåìåíè

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
Ôóíêöèÿ ðàñ÷¸òà áëèæàéøåãî ñîáûòèÿ äëÿ ÷àñòèöû

Àðãóìåíòû:
i - íîìåð ÷àñòèöû, äëÿ êîòîðîé ìû äîëæíû ðàññ÷èòàòü áëèæàéøåå ñîáûòèå
*/
void retime(long int &i) {
	particle p1 = particles[i];
	Box p1_box = boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]];
	long int jm = -1;  // ïåðåìåííàÿ äëÿ ñîõðàíåíèÿ òèïà áëèæàéøåãî ñîáûòèÿ
	double dt, dt_min = 1.0e+20;  // ïåðåìåííûå äëÿ ðàññ÷¸òà âðåìåíè äëÿ áëèæàéøåãî ñîáûòèÿ

	if (p1.vx < 0.0) {
		dt_min = (p1_box.x1 - p1.x[0]) / p1.vx;
		jm = -2;  // ñîáûòèå ïåðåñå÷åíèÿ ãðàíèöû Õ1 ÿ÷åéêè, â êîòîðîé íàõîäèòñÿ ÷àñòèöà 

				  // åñëè ìû íàõîäèìñÿ âáëèçè èäåàëüíîé ñòåíêè, òî ðàññ÷èòàòü âðåìÿ ñîóäàðåíèÿ
				  // ñ èäåàëüíîé ñòåíêîé
		if (p1.x_box[0] == 1) {
			dt_min = (p1_box.x1 + particle_R - p1.x[0]) / p1.vx;
			jm = -1;  // ñîáûòèå ñòîëêíîâåíèÿ ñ èäåàëüíîé ñòåíêîé
		}
	}
	else {
		dt_min = (p1_box.x2 - p1.x[0]) / p1.vx;
		jm = -4;  // ñîáûòèå ïåðåñå÷åíèÿ ãðàíèöû Õ2 ÿ÷åéêè, â êîòîðîé íàõîäèòñÿ ÷àñòèöà

		if (DIFFUSE_RIGHT_WALL == true) {
			dt = (p1.x_max - p1.x[0]) / p1.vx;
			if (dt < dt_min && dt > 0.0) {
				dt_min = dt;
				jm = -1;
			}

			dt = (L - particle_R - x_diffuse - p1.x[0]) / p1.vx;
			if (dt < dt_min && dt > 0.0) {
				dt_min = dt;
				jm = -20;
			}
		}
		else {
			// åñëè ìû íàõîäèìñÿ âáëèçè èäåàëüíîé ñòåíêè, òî ðàññ÷èòàòü âðåìÿ ñîóäàðåíèÿ
			// ñ èäåàëüíîé ñòåíêîé
			if (p1.x_box[0] == K2 - 1) {
				dt_min = (p1_box.x2 - particle_R - p1.x[0]) / p1.vx;
				jm = -1;  // ñîáûòèå ñòîëêíîâåíèÿ ñ èäåàëüíîé ñòåíêîé
			}
		}
	}

	if (p1.vy < 0.0) {
		dt = (p1_box.y1 - p1.y[0]) / p1.vy;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -5;  // ñîáûòèå ïåðåñå÷åíèÿ ãðàíèöû Y1 ÿ÷åéêè, â êîòîðîé íàõîäèòñÿ ÷àñòèöà  
		}
		if (p1.y_box[0] == 1) {
			dt = (p1_box.y1 + particle_Rx2 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -12;  // ñîáûòèå ðîæäåíèÿ îáðàçà ÷àñòèöû
			}
		}
		if (p1.y_box[0] == K - 1) {
			dt = (p1_box.y2 - particle_Rx2 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -11;  // ñîáûòèå óäàëåíèÿ îáðàçà ÷àñòèöû
			}
		}
	}
	else {
		dt = (p1_box.y2 - p1.y[0]) / p1.vy;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -6;  // ñîáûòèå ïåðåñå÷åíèÿ ãðàíèöû Y2 ÿ÷åéêè, â êîòîðîé íàõîäèòñÿ ÷àñòèöà
		}
		if (p1.y_box[0] == K - 1) {
			dt = (p1_box.y2 - particle_Rx2 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -11;  // ñîáûòèå ðîæäåíèÿ îáðàçà ÷àñòèöû
			}
		}
		if (p1.y_box[0] == 1) {
			dt = (p1_box.y1 + particle_Rx2 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -12;  // ñîáûòèå óäàëåíèÿ îáðàçà ÷àñòèöû
			}
		}
	}

	if (p1.vz < 0.0) {
		dt = (p1_box.z1 - p1.z[0]) / p1.vz;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -7;  // ñîáûòèå ïåðåñå÷åíèÿ ãðàíèöû Z1 ÿ÷åéêè, â êîòîðîé íàõîäèòñÿ ÷àñòèöà
		}
		if (p1.z_box[0] == 1) {
			dt = (p1_box.z1 + particle_Rx2 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -14;  // ñîáûòèå ðîæäåíèÿ îáðàçà ÷àñòèöû
			}
		}
		if (p1.z_box[0] == K - 1) {
			dt = (p1_box.z2 - particle_Rx2 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -13;  // ñîáûòèå óäàëåíèÿ îáðàçà ÷àñòèöû
			}
		}
	}
	else {
		dt = (p1_box.z2 - p1.z[0]) / p1.vz;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -8;  // ñîáûòèå ïåðåñå÷åíèÿ ãðàíèöû Z2 ÿ÷åéêè, â êîòîðîé íàõîäèòñÿ ÷àñòèöà
		}
		if (p1.z_box[0] == K - 1) {
			dt = (p1_box.z2 - particle_Rx2 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -13;  // ñîáûòèå ðîæäåíèÿ îáðàçà ÷àñòèöû
			}
		}
		if (p1.z_box[0] == 1) {
			dt = (p1_box.z1 + particle_Rx2 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -14;  // ñîáûòèå óäàëåíèÿ îáðàçà ÷àñòèöû
			}
		}
	}

	double temp, dx, dy, dz, dvx, dvy, dvz, d, dv, bij, dr;
	long int s, n, r, q, w, vp1 = 0, vp2 = -1;

	for (long int iv = 0; iv < 4; iv++) {
		if (p1.x[iv] > 0.0) {
			// Ïðîõîäèì ïî ÿ÷åéêàì, áëèæàéøèì ê ÿ÷åéêå, â êîòîðîé íàõîäèòñÿ ÷àñòèöà i
			for (r = p1.x_box[iv] - 1; r < p1.x_box[iv] + 2; r++)
				for (q = p1.y_box[iv] - 1; q < p1.y_box[iv] + 2; q++)
					for (w = p1.z_box[iv] - 1; w < p1.z_box[iv] + 2; w++) {

						// åñëè èíäåêñ ÿ÷åéêè âûõîäèò çà ãðàíèöó ñèñòåìû òî
						// ïåðåõîäèì íà ñëåäóþùèé øàã öèêëà
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

						// Ïðîõîäèì ïî âñåì ÷àñòèöàì â âûáðàííîé ÿ÷åéêå è ïðîâåðÿåì âîçìîæíîñòü
						// ñòîëêíîâåíèÿ ýòèõ ÷àñòèö ñ ÷àñòèöåé i
						for (s = 0; s <= boxes_yz[q][w][r].end; ++s) {
							n = boxes_yz[q][w][r].particles[s];
							long int jv = 0;

							if (n >= NP) {
								jv = n / NP; // íîìåð èçíà÷àëüíîé ÷àñòèöû
								n = n % NP;  // íîìåð îáðàçà
							}

							/*
							Cîõðàíÿåì â ïåðåìåííóþ p âñå äàííûå ÷àñòèöû n,
							÷òîáû äàëåå çàïèñûâàòü âñå îïåðàöèè êîðî÷å è áåç îáðàùåíèÿ
							â ãëîáàëüíóþ ïàìÿòü
							*/
							particle p = particles[n];

							/*
							Ðàññ÷èòûâàåì ðàçíèöó ñîáñòâåííîãî âðåìåíè äâóõ ÷àñòèö,
							ýòî íåîáõîäèìî ÷òîáû ñèíõðîíèçîâàòü èõ ìåæäó ñîáîé è ðàññ÷èòûâàòü
							âîçìîæíîñòü è âðåìÿ ñòîëêíîâåíèÿ â ñèñòåìå, ãäå îáå ÷àñòèöû áóäóò
							èìåòü îäèíàêîâîå ñîáñòâåííîå âðåìÿ
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
								dr = particle_D2 - dx * dx - dy * dy - dz * dz;
								// ðàññ÷èòûâàåì äèñêðèìèíàíò â óðàâíåíèè äëÿ âû÷èñëåíèÿ âðåìåíè ñîóäàðåíèÿ
								d = bij * bij + dv * dr;

								// åñëè äèñêðèìèíàíò áîëüøå íóëÿ òî ñîóäàðåíèå âîçìîæíî
								if (d > 0.0) {
									dt = -(sqrt(d) + bij) / dv;

									/*
									Ê ðàçíèöå â ñîáñòâåííîì âðåìåíè ÷àñòèö ïðèáàâëÿåì âðåìÿ äî èõ ñîóäàðåíèÿ,
									â ðåçóëüòàòå ïîëó÷àåì âðåìÿ dt ÷åðåç êîòîðîå ýòî ñîóäàðåíèå ïðîèçîéä¸ò äëÿ
									÷àñòèöû p, òàêèì îáðàçîì ìû ïîëó÷àåì âîçìîæíîñòü ñðàâíèòü âðåìÿ òîëüêî ÷òî
									ðàññ÷èòàííîãî ñîóäàðåíèÿ ñî âðåìåíåì áëèæàéøåãî ñîáûòèÿ äëÿ ÷àñòèöû p è
									ïðèíÿòü ðåøåíèå - äîëæíû ëè ìû ïåðåçàïèñàòü ñîáûòèå äëÿ ÷àñòèöû p èëè íåò.
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
	Åñëè íîâîå ñîáûòèå - ýòî ñòîëêíîâåíèå äâóõ ÷àñòèö, òî äëÿ âòîðîé ÷àñòèöû
	íåîáõîäèìî ïåðåçàïèñàòü å¸ ïðåäûäóùåå áëèæàéøåå ñîáûòèå è ñèíõðîíèçîâàòü
	ýòó ÷àñòèöó ïî âðåìåíè ñ ÷àñòèöåé i
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
	if (i == 168 || jm == 168 || i == 171 || jm == 171) {
	printf("\n retime result: %d %d %d %d %.15le", i, jm, vp1, vp2, dt_min);
	printf("\n 168[0] x_box, y_box, z_box %d %d %d", particles[168].x_box[0], particles[168].y_box[0], particles[168].z_box[0]);
	printf("\n 171[3] x_box, y_box, z_box %d %d %d \n", particles[171].x_box[3], particles[171].y_box[3], particles[171].z_box[3]);
	}
	*/

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
	Â ñëó÷àå åñëè ïðè ðàñ÷¸òå áëèæàéøåãî ñîáûòèÿ ìû ïîëó÷èëè îòðèöàòåëüíîå âðåìÿ
	íåîáõîäèìî ïðåðâàòü âûïîëíåíèå ïðîãðàììû è ðàñïå÷àòàòü îòëàäî÷íóþ èíôîðìàöèþ
	î ðàññ÷èòàííîì ñîáûòèè
	*/
	if (dt_min < -1.0e-11) {
		printf("\n retime result: %ld %ld, %.16le\n ", i, jm, dt_min);
		printf("\n p1.x = %.16le, p1.y = %.16le, p1.z = %.16le \n", p1.x[0], p1.y[0], p1.z[0]);
		printf("\n p1.vx = %.16le, p1.vy = %.16le, p1.vz = %.16le \n", p1.vx, p1.vy, p1.vz);
		printf("\n p1.x_box = %ld, p1.y_box = %ld, p1.z_box = %ld", p1.x_box[0], p1.z_box[0], p1.y_box[0]);
		printf("\n im = %ld, jm = %ld, dt = %.16le, A = %.16le", i, jm, dt_min, A);
		printf("\n p1.box.x = [%.16le; %.16le]", boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].x1, boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].x2);
		printf("\n p1.box.y = [%.16le; %.16le]", boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].y1, boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].y2);
		printf("\n p1.box.z = [%.16le; %.16le]", boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].z1, boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]].z2);
		throw "dt < 0!";
	}
}


void change_particles(long int i, long int vp1, long int vp2) {
	Box box1, box2;
	long int y_box1, z_box1, x_box1, y_box2, z_box2;
	long int j;

	/*
	if (particles[i].x[vp1] < 0.0) {
	printf("AAAAAAAA");
	throw "A";
	}
	if (particles[i].x[vp2] < 0.0) {
	printf("BBBBBBB");
	throw "B";
	}
	*/

	x_box1 = particles[i].x_box[vp1];
	y_box1 = particles[i].y_box[vp1];
	z_box1 = particles[i].z_box[vp1];

	box1 = boxes_yz[y_box1][z_box1][x_box1];
	box1.particles[particles[i].box_i[vp1]] = i + NP * vp2;
	boxes_yz[y_box1][z_box1][x_box1] = box1;

	y_box2 = particles[i].y_box[vp2];
	z_box2 = particles[i].z_box[vp2];

	box2 = boxes_yz[y_box2][z_box2][x_box1];
	box2.particles[particles[i].box_i[vp2]] = i + NP * vp1;
	boxes_yz[y_box2][z_box2][x_box1] = box2;

	j = particles[i].box_i[vp1];
	particles[i].box_i[vp1] = particles[i].box_i[vp2];
	particles[i].box_i[vp2] = j;

	particles[i].y_box[vp1] = y_box2;
	particles[i].z_box[vp1] = z_box2;
	particles[i].y_box[vp2] = y_box1;
	particles[i].z_box[vp2] = z_box1;
}

/*
Ôóíêöèÿ îáìåíà ÷àñòèöû è "îáðàçà" ïðè ïåðåñå÷åíèè ÷àñòèöåé ïåðèîäè÷åñêèõ
ãðàíè÷íûõ óñëîâèé.

Àðãóìåíòû:
im - íîìåð ÷àñòèöû, ïåðåñåêàþùåé ïåðèîäè÷åñêèå ãðàíè÷íûå óñëîâèÿ
jm - íîìåð ãðàíèöû, ÷åðåç êîòîðóþ ïðîõîäèò öåíòð ÷àñòèöû
*/
void change_with_virt_particles(long int &im, long int &jm) {
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
			k = particles[im].z[1];
			particles[im].z[1] = particles[im].z[3];
			particles[im].z[3] = k;

			change_particles(im, 1, 3);
		}
	}

}


/*
Ôóíêöèÿ ñîçäàíèÿ íîâîãî "îáðàçà".

Àðãóìåíòû:
i - íîìåð ÷àñòèöû, îáðàç êîòîðîé íåîáõîäèìî ñîçäàòü
need_to_check - ôëàã, ïîçâîëÿþùèé íå ïðîâåðÿòü íóæåí ëè îáðàç äëÿ äàííîé
÷àñòèöû èëè íåò, è ñðàçó ñîçäàâàòü îáðàç (êîãäà ìû óâåðåíû,
÷òî îáðàç íåîáõîäèìî ñîçäàòü).
*/
void create_virt_particle(long int &im, long int &jm) {
	Box box1;

	switch (jm)
	{
	case -11:     // ÷àñòèöà íà ãðàíèöå Y = À-2R
		particles[im].y[1] = -particle_Rx2;
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
	case -12:     // ÷àñòèöà íà ãðàíèöå Y = 2R
		particles[im].y[1] = A + particle_Rx2;
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
	case -13:     // ÷àñòèöà íà ãðàíèöå Z = A-2R
		particles[im].z[2] = -particle_Rx2;
		particles[im].x[2] = particles[im].x[0];
		particles[im].y[2] = particles[im].y[0];
		particles[im].x_box[2] = particles[im].x_box[0];
		particles[im].z_box[2] = 0;
		particles[im].y_box[2] = particles[im].y_box[0];

		box1 = boxes_yz[particles[im].y_box[2]][0][particles[im].x_box[2]];
		box1.end++;
		box1.particles[box1.end] = im + NP * 2;
		boxes_yz[particles[im].y_box[2]][0][particles[im].x_box[2]] = box1;
		particles[im].box_i[2] = box1.end;

		break;
	case -14:     // ÷àñòèöà íà ãðàíèöå Z = 2R
		particles[im].z[2] = A + particle_Rx2;
		particles[im].x[2] = particles[im].x[0];
		particles[im].y[2] = particles[im].y[0];
		particles[im].x_box[2] = particles[im].x_box[0];
		particles[im].z_box[2] = K;
		particles[im].y_box[2] = particles[im].y_box[0];

		box1 = boxes_yz[particles[im].y_box[2]][K][particles[im].x_box[2]];
		box1.end++;
		box1.particles[box1.end] = im + NP * 2;
		boxes_yz[particles[im].y_box[2]][K][particles[im].x_box[2]] = box1;
		particles[im].box_i[2] = box1.end;

		break;
	}

	// åñëè ó ÷àñòèöû åñòü äâà îáðàçà, òî íåîáõîäèìî ñîçäàòü òðåòèé îáðàç
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
		box1.particles[box1.end] = im + NP * 3;
		boxes_yz[y][z][particles[im].x_box[3]] = box1;
		particles[im].box_i[3] = box1.end;
	}
}


/*
Ôóíêöèÿ çàãðóçêè èíôîðìàöèè î íåêîòîðîé ÷àñòèöå èç ôàéëà.
Ïðè çàãðóçêå ñèñòåìû ìû ïîñëåäîâàòåëüíî ñ÷èòûâàåì èíôîðìàöèþ
î êàæäîé ÷àñòèöå ñ ïîìîùüþ ýòîé ôóíêöèè.

Àðãóìåíòû:
loading_file - ññûëêà íà ôàéë, èç êîòîðîãî íåîáõîäèìî ñ÷èòàòü
äàííûå î ÷àñòèöå èëè îáðàçå.
*/
void load_information_about_one_particle(FILE *loading_file) {
	double a1, a2, a3;
	long i, x_box, y_box, z_box, end;

	// Ñ÷èòûâåì íîìåð ÷àñòèöû
	fscanf(loading_file, "%ld", &i);

	// Ñ÷èòûâàåì êîîðäèíàòû ÷àñòèöû
	fscanf(loading_file, "%le %le %le", &a1, &a2, &a3);
	particles[i].x[0] = a1;
	particles[i].y[0] = a2;
	particles[i].z[0] = a3;

	// Ñ÷èòûâàåì èíôîðìàöèþ î ñêîðîñòÿõ ÷àñòèöû
	fscanf(loading_file, "%le %le %le", &a1, &a2, &a3);
	particles[i].vx = a1;
	particles[i].vy = a2;
	particles[i].vz = a3;

	particles[i].t = 0.0;
	particles[i].dt = 0.0;
	particles[i].ti = 1;

	if (particles[i].y[0] < particle_Rx2) {
		particles[i].x[1] = particles[i].x[0];
		particles[i].y[1] = particles[i].y[0] + A;
		particles[i].z[1] = particles[i].z[0];
	}
	if (particles[i].y[0] > A - particle_Rx2) {
		particles[i].x[1] = particles[i].x[0];
		particles[i].y[1] = particles[i].y[0] - A;
		particles[i].z[1] = particles[i].z[0];
	}
	if (particles[i].z[0] < particle_Rx2) {
		particles[i].x[2] = particles[i].x[0];
		particles[i].z[2] = particles[i].z[0] + A;
		particles[i].y[2] = particles[i].y[0];
	}
	if (particles[i].z[0] > A - particle_Rx2) {
		particles[i].x[2] = particles[i].x[0];
		particles[i].z[2] = particles[i].z[0] - A;
		particles[i].y[2] = particles[i].y[0];
	}
	if (particles[i].x[1] > 0.0 && particles[i].x[2] > 0.0) {
		particles[i].x[3] = particles[i].x[0];
		particles[i].z[3] = particles[i].z[2];
		particles[i].y[3] = particles[i].y[1];
	}

	for (long int j = 0; j < 4; j++) {
		if (particles[i].x[j] > 0.0) {
			// Îïðåäåëÿåì â êàêóþ ÿ÷åéêó íåîáõîäèìî çàïèñàòü íîâóþ ÷àñòèöó
			x_box = long((particles[i].x[j] + dL) / dL);
			y_box = long((particles[i].y[j] + dA) / dA);
			z_box = long((particles[i].z[j] + dA) / dA);

			if (particles[i].x[j] < boxes_yz[y_box][z_box][x_box].x1 ||
				particles[i].x[j] > boxes_yz[y_box][z_box][x_box].x2) {
				printf("\n %ld \n", long int((particles[i].x[j] + dL) / dL));
				printf("%d %d %d", x_box, y_box, z_box);
				throw "aaa";
			}


			if (y_box < 0) y_box = 0;
			if (z_box < 0) z_box = 0;
			if (y_box > K) y_box = K;
			if (z_box > K) z_box = K;

			particles[i].x_box[j] = x_box;
			particles[i].y_box[j] = y_box;
			particles[i].z_box[j] = z_box;

			// çàïèñûâàåì ÷àñòèöó â ÿ÷åéêó
			end = particles[i].box_i[j] = ++boxes_yz[y_box][z_box][x_box].end;
			boxes_yz[y_box][z_box][x_box].particles[end] = i + NP*j;

			/*
			Åñëè â îäíîé ÿ÷åéêå íàõîäèòñÿ ñëèøêîì ìíîãî ÷àñòèö,
			òî çàâåðøèòü ïðîãðàììó, âûâåäÿ ñîîáùåíèå îá îøèáêå.
			*/
			if (end > 10) {
				printf("Error: Too many particles in one box");
				printf("\n x_box, y_box, z_box = %ld %ld %ld , end = %ld \n", x_box, y_box, z_box, end);
				throw "Error: Too many particles in one box";
			}
		}
	}

	global_E += a1*a1 + a2*a2 + a3*a3;

	/*
	Ïðîâåðÿåì ÷òî â äàííîé ÿ÷åéêå íåò ïîâòîðÿþùèõñÿ íîìåðîâ ÷àñòèö
	*/
	Box b = boxes_yz[particles[i].y_box[0]][particles[i].z_box[0]][particles[i].x_box[0]];
	for (long int t = 0; t < particles[i].box_i[0]; ++t)
		for (long int d = t + 1; d <= particles[i].box_i[0]; ++d) {
			if (b.particles[t] == b.particles[d])
			{
				printf("\n Duplicated particles in one box: %ld \n", i);

				for (long int j = 0; j < b.end; ++j)
					printf("%ld ", b.particles[j]);

				throw "Duplicated particles in one box";
			}
		}
}


/*
Ôóíêöèÿ çàãðóçêè äàííûõ î ñèñòåìå èç ôàéëà.
Çàãðóæàþòñÿ ïàðàìåòðû îáú¸ìà, êîîðäèíàòû è ñêîðîñòè âñåõ ÷àñòèö.

Àðãóìåíòû:
file_name - èìÿ ôàéëà, ñîäåðæàùåãî âñå äàííûå î ñèñòåìå.
*/
void load_seed(std::string file_name) {
	double a1, x, y, z;
	long h;
	FILE *loading_file;

	global_E = 0.0;

	/*
	Îòêðûâàåì óêàçàííûé ôàéë äëÿ ÷òåíèÿ,
	ïèøåì íà ýêðàí èíôîðìàöèþ îá îøèáêå è âûõîäèì èç ïðîãðàììû åñëè
	óêàçàííîãî ôàéëà íå ñóùåñòâóåò èëè åãî íåâîçìîæíî îòêðûòü äëÿ ÷òåíèÿ.
	*/
	try {
		loading_file = fopen(file_name.c_str(), "r");
	}
	catch (...) {
		printf("\n ERROR: Can't read file '%s' - is it exist? \n", file_name.c_str());
		exit(1);
	}

	// ñ÷èòûâàåì ïëîòíîñòü
	fscanf(loading_file, "%le", &a1);

	// ñ÷èòûâàåì êîëè÷åñòâî ÷àñòèö â ñèñòåìå
	fscanf(loading_file, "%ld", &h);
	NP = h;
	// ñ÷èòûâàåì çíà÷åíèå A
	fscanf(loading_file, "%le", &a1);
	A = a1;
	// ñ÷èòûâàåì çíà÷åíèå L 
	fscanf(loading_file, "%le", &a1);
	L = a1;
	// ñ÷èòûâàåì çíà÷åíèå R 
	//fscanf(loading_file, "%le", &a1);
	//change_R(a1);

	// EN: search for the appropriate values for dA, dL, K, K2
	// RU: èùåì ïîäõîäÿùåå çíà÷åíèå äëÿ dA, dL, K, K2
	dA = 2.2*particle_R;
	dL = 2.2*particle_R;
	K = long(A / dA) + 1;  // êîëè÷åñòâî ÿ÷ååê ïî Y è Z
	K2 = long(L / dL) + 1; // êîëè÷åñòâî ÿ÷ååê ïî X
	dA = A / (K - 1);       // òî÷íûå ðàçìåðû ÿ÷ååê ïî X, Y è Z
	dL = L / (K2 - 1);      // âûáðàííûå òàê, ÷òîáû, ÿ÷åéêè ñîâïàäàëè
							// ñ ðàçìåðàìè îáú¸ìà

							/*
							Ïîñëå òîãî êàê ìû ðàññ÷èòàëè êîëè÷åñòâî ÿ÷ååê è èõ ðàçìåðû
							íåîáõîäèìî èíèöèàëèçèðîâàòü äàííûå êàæäîé ÿ÷åéêè, óêàçàòü
							ãðàíèöû êàæäîé ÿ÷åéêè ïî X, Y, Z, ò.ê. ýòè äàííûå áóäóò
							èñïîëüçîâàòüñÿ äëÿ ðàññ÷¸òà âðåìåíè äî ïåðåñå÷åíèÿ ãðàíèö
							ÿ÷ååê ÷àñòèöàìè, êîòîðûå íàõîäÿòñÿ â ýòèõ ÿ÷åéêàõ.
							*/
	y = z = -dA;
	x = -dL;
	for (long i = 0; i <= K; i++) {
		for (long j = 0; j <= K; j++) {
			for (long w = 0; w <= K2; w++) {
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
	Îáíóëÿåì âñå äàííûå äëÿ âñåõ ÷àñòèö è îáðàçîâ ïåðåä ñ÷èòûâàíèåì
	ôàéëà ñ äàííûìè.
	Äëÿ îáðàçîâ, êîòîðûå íå îïèñàíû â ôàéëå ñîõðàíåíèÿ, îñòàíóòñÿ ýòè
	íà÷àëüíûå çíà÷åíèÿ.
	*/
	for (long i = 0; i < NP; i++) {

		for (int j = 0; j < 4; j++) {
			particles[i].x[j] = -10.0*particle_R;
			particles[i].y[j] = -10.0*particle_R;
			particles[i].z[j] = -10.0*particle_R;
			particles[i].x_box[j] = 0;
			particles[i].y_box[j] = 0;
			particles[i].z_box[j] = 0;
			particles[i].box_i[j] = 0;
		}

		particles[i].ti = 0;
		particles[i].dt = 0.0;
	}

	/*
	Çàãðóæàåì äàííûå î âñåõ ÷àñòèöàõ è "îáðàçàõ" èç ôàéëà ñîõðàíåíèÿ
	*/
	for (long i = 0; i < NP; i++) {
		load_information_about_one_particle(loading_file);

		if (DIFFUSE_RIGHT_WALL == true && particles[i].x[0] >= L - particle_R - x_diffuse) {
			if (particles[i].vx < 0.0) {
				particles[i].x_max = particles[i].x[0];
			}
			else {
				double dx = (L - particle_R - particles[i].x[0]) / (rand() / 100.0 + 10.0);
				particles[i].x_max = particles[i].x[0] + dx;

				if (dx > x_diffuse || dx < 0.0 || particles[i].x_max < particles[i].x[0] || particles[i].x_max > L - particle_R) {
					printf("\n dx = %.15le \n", dx);
					throw "alarm";
				}
			}
		}
		else {
			particles[i].x_max = L - particle_R;
		}
	}

	fclose(loading_file);

	printf("DF");

	/*
	Èíèöèàëèçèðóåì î÷åðåäü ñîáûòèé ïóñòûìè ñîáûòèÿìè è óñòàíàâëèâàåì
	óêàçàòåëü êîíöà ñïèñêà íà åãî íà÷àëî.
	*/
	last = 1;
	time_queue[0].t = 0.0;
	time_queue[0].im = -1;
	for (long i = 1; i < NP; i++) time_queue[i].t = 1.0E+20;

	/*
	Ðàññ÷èòûâàåì íîâûå ñîáûòèÿ äëÿ âñåõ ÷àñòèö è ñóùåñòâóþùèõ â ñèñòåìå "îáðàçîâ"
	*/
	for (long i = 0; i < NP; i++) {
		retime(i);
	}
}


/*
Ôóíêöèÿ ñîõðàíåíèÿ ñîñòîÿíèÿ ñèñòåìû â òåêñòîâûé ôàéë.
ïåðåä ñîõðàíåíèåì ïðîèçâîäèòñÿ ñèíõðîíèçàöèÿ ÷àñòèö ïî âðåìåíè.

Àðãóìåíòû:
file_name - èìÿ ôàéëà, â êîòîðûé áóäåò çàïèñàíà èíôîðìàöèÿ
*/
void save(std::string file_name) {
	double x, y, z, dt;
	double t_global = get_maximum_particle_time();

	FILE *save_file = fopen(file_name.c_str(), "w+");

	double etta = (4.0 * PI * particle_R3 * NP) / (3.0 * A * A * (L - particle_Rx2));
	fprintf(save_file, "%.15le\n", etta);

	// RU: ñîõðàíÿåì èíôîðìàöèþ î êîëè÷åñòâå ÷àñòèö è ðàçìåðàõ ñèñòåìû
	// EN: save information about count of particles and size of the system
	fprintf(save_file, "%ld\n", NP);
	fprintf(save_file, "%.15le\n", A);
	fprintf(save_file, "%.15le\n", L);
	//fprintf(save_file, "%.15le\n", particle_R);

	// RU: ñîõðàíÿåì êîîðäèíàòû âñåõ ÷àñòèö è èõ ñêîðîñòè: x, y, z, vx, vy, vz
	// EN: we need to save all coordinates of particles: x, y, z, vx, vy, vz
	for (long i = 0; i < NP; ++i) {
		fprintf(save_file, "%ld\n", i);

		/*
		Ñèíõðîíèçèðóåì ÷àñòèöó i ñ ãëîáàëüíûì âðåìåíåì ñèñòåìû,
		â èòîãå äàííûå âñåõ ÷àñòèö áóäóò çàïèñàíû äëÿ îäíîãî ìîìåíòà âðåìåíè,
		êîòîðûé ñîâïàäàåò ñî âðåìåíåì ïîñëåäíåãî ïðîèçîøåäøåãî â ñèñòåìå ñîáûòèÿ.
		*/
		dt = t_global - particles[i].t;

		/*
		Ïðè ñîõðàíåíèè ñîñòîÿíèÿ ñèñòåìû ìû óêàçûâàåì
		àáñîëþòíûå êîîðäèíàòû ÷àñòèö è ïàðàìåòðû ñèñòåìû, òàê, ÷òîáû
		öåíòð ñèñòåìû áûë â òî÷êå (L/2, A/2), à íå â (0, 0).
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
Ïðîöåäóðà íîâîãî "ïîñåâà" ñèñòåìû.

Àðãóìåíòû:
NN - ÷èñëî ÷àñòèö â ðåáðå îáú¸ìî öåíòðèðîâàííîãî êðèñòàëëà, íà îñíîâå
êîòîðîãî äåëàåòñÿ èçíà÷àëüíûé ïîñåâ
etta - íà÷àëüíàÿ ïëîòíîñòü, êîòîðóþ íåîáõîäèìî çàäàòü â ñèñòåìå
lamda - ðàññòîÿíèå ìåæäó óçëàìè
*/
void new_seed(long particles_count, double etta) {
	long i = 0, j, j1, j2, j3;

	// default 16x48
	double lamda = 2.0*particle_R + 0.1*particle_R;
	long Qyz = particles_count;  // ÷èñëî óçëîâ ïî ó è z
	long Qx = particles_count * 4; // ÷èñëî óçëîâ ïî x
	double vx, vy, vz;
	double delta_L;

	NP = Qx * Qyz * Qyz;
	A = Qyz * lamda;
	L = ((4.0 * PI * particle_R3 * NP) / etta) / (3.0 * A * A) + 2.0*particle_R;
	delta_L = L / Qx;

	double x = delta_L / 2.0;
	double y = lamda / 2.0;
	double z = lamda / 2.0;

	if (x < particle_R || y < particle_R) {
		printf("\n Incorrect parameters for new seed! \n");
		if (x < particle_R) {
			printf("Delta X is to small = %.5le, should be > %.5le. Please change etta0 \n", x, particle_R);
		}
		else {
			printf("Delta Y is to small = %.5le, should be > %.5le. Please change lamda \n", y, particle_R);
		}
		throw "incorrect parameters";
	}

	//srand(time(NULL));  - ñîâñåì ñëó÷àéíî
	srand(10);   // - íå ñëó÷àéíî

	while (i < NP) {

		vx = double(double(rand() % 30000 + 100) / (30000.0 + double(rand() % 1000))
			- double(rand() % 30000 + 100) / (30000.0 + double(rand() % 1000)));
		vy = double(double(rand() % 30000 + 100) / (30000.0 + double(rand() % 1000))
			- double(rand() % 30000 + 100) / (30000.0 + double(rand() % 1000)));
		vz = double(double(rand() % 30000 + 100) / (30000.0 + double(rand() % 1000))
			- double(rand() % 30000 + 100) / (30000.0 + double(rand() % 1000)));

		for (j1 = 0; j1 < 2; j1++) {
			for (j2 = 0; j2 < 2; j2++) {
				for (j3 = 0; j3 < 2; j3++) {
					j = i + j1 * 4 + j2 * 2 + j3;

					particles[j].x[0] = j1*L + x - 2.0*x*j1; // x or L-x
					particles[j].y[0] = j2*A + y - 2.0*y*j2; // y or A-y
					particles[j].z[0] = j3*A + z - 2.0*z*j3; // z or A-z

					if ((j1 + j2 + j3) % 2 == 0) {
						particles[j].vx = -vx;
						particles[j].vy = vy;
						particles[j].vz = -vz;
					}
					else {
						particles[j].vx = vx;
						particles[j].vy = -vy;
						particles[j].vz = vz;
					}

					particles[j].t = 0.0;
					particles[j].dt = 0.0;
					particles[j].ti = 0;
				}
			}
		}
		i += 8;
		y += lamda;
		if (y > A / 2) {
			y = lamda / 2;
			z += lamda;
			if (z > A / 2) {
				z = lamda / 2;
				x += delta_L;
			}
		}
	}

	double Lx = 0.0;
	for (i = 0; i < NP; i++) {
		Lx += (A / 2.0 - particles[i].y[0]) * particles[i].vz - (A / 2.0 - particles[i].z[0]) * particles[i].vy;
	}
	printf("\n Lx = %.15le\n", Lx);

	/*
	Cîõðàíÿåì è çàãðóæàåì ñèñòåìó èç ôàéëà ÷òîáû èíèöèàëèçèðîâàòü
	âñå íåîáõîäèìûå ïåðåìåííûå.
	*/
	save("new.txt");
	load_seed("new.txt");
}


void clear_cell_from_virt_particle(long i, long j) {
	particle p1 = particles[i];
	Box p1_box = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]];

	if (p1_box.particles[p1.box_i[j]] % NP == i) {
		long end = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].end;

		if (p1.box_i[j] <= end) {
			/*
			Ñòèðàåì ÷àñòèöó èç ÿ÷åéêè â êîòîðîé îíà íàõîäèëàñü,
			âìåñòî íå¸ âñòàâëÿåì ÷àñòèöó èç êîíöà ñïèñêà ÷àñòèö â ÿ÷åéêå è
			óìåíüøàåì ÷èñëî ÷àñòèö â ÿ÷åéêå íà 1.
			*/
			long fj = p1_box.particles[end] % NP;
			long kj = p1_box.particles[end] / NP;

			p1_box.particles[p1.box_i[j]] = p1_box.particles[end];
			particles[fj].box_i[kj] = p1.box_i[j];
			--p1_box.end;
			boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]] = p1_box;
		}
	}
}


/*
Ôóíêöèÿ èçìåíåíèÿ ñîñòîÿíèÿ ÷àñòèö â ñîîòâåòñòâèè ñ íàñòóïèâøèì
â ñèñòåìå ñîáûòèåì.

Àðãóìåíòû:
im - íîìåð ÷àñòèöû
jm - íîìåð âòîðîé ÷àñòèöû èëè íîìåð ãðàíèöû
*/
void reform(long &im, long &jm, long &vp1, long &vp2) {
	particle p1 = particles[im];
	double q1, q2, z, dx, dy, dz;

	/*
	Åñëè jm >= 0 òî íàñòóïèâøåå ñîáûòèå - ýòî ñîóäàðåíèå äâóõ
	÷àñòèö èëè îáðàçîâ.
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

		q1 = (dx * p1.vx + dy * p1.vy + dz * p1.vz) / particle_D2;
		q2 = (dx * p2.vx + dy * p2.vy + dz * p2.vz) / particle_D2;
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
		switch (jm) {
		case -1:   // ñîóäàðåíèå ñ èäåàëüíîé ñòåíêîé
			particles[im].vx = -particles[im].vx;

			if (particles[im].x[0] < L / 2.0) {
				dV_left_wall += abs(particles[im].vx);
				N_collissions_left_wall += 1;
			}
			else {
				dV_right_wall += abs(particles[im].vx);
				N_collissions_right_wall += 1;
			}

			break;
		case -11:  // ïåðåñå÷åíèå ãðàíèöû Y = A - 2R
			if (p1.vy > 0.0) {
				create_virt_particle(im, jm);
			}
			else {
				particles[im].x[1] = -10.0*particle_R;
				particles[im].x[3] = -10.0*particle_R;

				clear_cell_from_virt_particle(im, 1);
				clear_cell_from_virt_particle(im, 3);
			}
			break;
		case -12:  // ïåðåñå÷åíèå ãðàíèöû Y = 2R
			if (p1.vy < 0.0) {
				create_virt_particle(im, jm);
			}
			else {
				particles[im].x[1] = -10.0*particle_R;
				particles[im].x[3] = -10.0*particle_R;

				clear_cell_from_virt_particle(im, 1);
				clear_cell_from_virt_particle(im, 3);
			}
			break;
		case -13:    // ïåðåñå÷åíèå ãðàíèöû Z = A - 2R
			if (p1.vz > 0.0) {
				create_virt_particle(im, jm);
			}
			else {
				particles[im].x[2] = -10.0*particle_R;
				particles[im].x[3] = -10.0*particle_R;

				clear_cell_from_virt_particle(im, 2);
				clear_cell_from_virt_particle(im, 3);
			}
			break;
		case -14:  // ïåðåñå÷åíèå ãðàíèöû Z = 2R
			if (p1.vz < 0.0) {
				create_virt_particle(im, jm);
			}
			else {
				particles[im].x[2] = -10.0*particle_R;
				particles[im].x[3] = -10.0*particle_R;

				clear_cell_from_virt_particle(im, 2);
				clear_cell_from_virt_particle(im, 3);
			}
			break;
		case -20:
			particles[im].x_max = L - particle_R - 0.999 * x_diffuse + (0.998*x_diffuse*0.0001)*(rand() % 10000);
			break;
		case -100:
			break;
		default:
			for (int j = 0; j < 4; ++j) {
				if (p1.x[j] > 0.0) {
					Box p1_box = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]];
					long int end = p1_box.end;

					/*
					Ñòèðàåì ÷àñòèöó èç ÿ÷åéêè â êîòîðîé îíà íàõîäèëàñü,
					âìåñòî íå¸ âñòàâëÿåì ÷àñòèöó èç êîíöà ñïèñêà ÷àñòèö â ÿ÷åéêå è
					óìåíüøàåì ÷èñëî ÷àñòèö â ÿ÷åéêå íà 1.
					*/
					long int fj = p1_box.particles[end] % NP;
					long int kj = p1_box.particles[end] / NP;

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

					switch (jm) {
					case -2:
						--p1.x_box[j];
						if (p1.x[j] > 0.0) {
							p1.x[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].x2;
						}
						break;
					case -4:
						++p1.x_box[j];
						if (p1.x[j] > 0.0) {
							p1.x[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].x1;
						}
						break;
					case -5:
						--p1.y_box[j];
						p1.y[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].y2;
						break;
					case -6:
						++p1.y_box[j];
						p1.y[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].y1;
						break;
					case -7:
						--p1.z_box[j];
						p1.z[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].z2;
						break;
					case -8:
						++p1.z_box[j];
						p1.z[j] = boxes_yz[p1.y_box[j]][p1.z_box[j]][p1.x_box[j]].z1;
						break;
					}

					/*
					Ïîñëå òîãî, êàê èíäåêñ ÿ÷åéêè äëÿ ÷àñòèöû áûë èçìåí¸í,
					äîáàâëÿåì ÷àñòèöó â íîâóþ ÿ÷åéêó (èëè âîçâðàùàåì å¸ â ñòàðóþ
					ÿ÷åéêó, åñëè èíäåêñ ÿ÷åéêè íå èçìåíèëñÿ).
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
			break;
		}
	}
}


void check_overlap() {
	double time = get_maximum_particle_time();

	for (long int i = 0; i < particles_for_check_count; i++) {
		for (long int j = i + 1; j < particles_for_check_count; j++) {
			particle p01 = particles[particles_for_check[i]];
			particle p02 = particles[particles_for_check[j]];

			double dt01 = time - p01.t;
			double dt02 = time - p02.t;
			double ddx = p01.x[0] + p01.vx*dt01 - p02.x[0] - p02.vx*dt02;
			double ddy = p01.y[0] + p01.vy*dt01 - p02.y[0] - p02.vy*dt02;
			double ddz = p01.z[0] + p01.vz*dt01 - p02.z[0] - p02.vz*dt02;
			double rr = ddx*ddx + ddy*ddy + ddz*ddz;

			if (rr < particle_D2 - 1.0e-14) {
				printf("\n r = %.15le ", rr);
				printf("\n dt01 = %.15le < p01.dt = %.15le ? \n ", dt01, p01.dt);
				printf("\n dt02 = %.15le < p02.dt = %.15le ? \n ", dt02, p02.dt);

				throw "ALARM";
			}
		}
	}
}


/*
Ôóíêöèÿ "øàã", îñíîâíîé öèêë ïðîãðàììû,
ïðîèçâîäèò ðàññ÷¸ò äèíàìèêè ñèñòåìû â òå÷åíèè 1 ñîóäàðåíèÿ
íà êàæäóþ ÷àñòèöó  - (NP/2) ñòîëêíîâåíèé â ñèñòåìå.
*/
void step() {
	particle p1;
	long int i, im, jm, vp1, vp2;
	double time = get_maximum_particle_time();

	COLL_COUNT = 0;
	jm = 0;

	while (COLL_COUNT < NP / 2 || jm < 0) {
		// ñ÷èòûâàåì ïåðâîå ñîáûòèå èç ëèíåéêè ñîáûòèé
		im = time_queue[1].im;
		jm = time_queue[1].jm;
		vp1 = time_queue[1].vp1;
		vp2 = time_queue[1].vp2;

		// óäàëÿåì ïåðâîå ñîáûòèå
		delete_event(1);

		p1 = particles[im];

		/*
		Îáíóëÿåì óêàçàòåëü íà ñîáûòèå ÷àñòèö, ó÷àâñòâóþùèõ â ýòîì
		ñîáûòèè, ïðèðàâíèâàÿ èõ -1, òàê, ÷òîáû îíè íå óêàçûâàëè íà
		ñóùåñòâóþùèå ñîáûòèÿ.
		*/
		p1.ti = -1;
		if (jm >= 0) particles[jm].ti = -1;

		/*
		Åñëè íàñòóïèëî ñîáûòèå -100, çíà÷èò ïðîñòî
		ñèíõðîíèçèðóåì ÷àñòèöó ñ êîòîðîé ïðîèçîøëî
		äàííîå ñîáûòèå ñ ãëîáàëüíûì âðåìåíåì ñèñòåìû.
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

		dP_time += p1.dt;

		p1.dt = 0.0;
		time = p1.t;


		particles[im] = p1;

		// Ïðîèçâîäèì èçìåíåíèÿ â ñèñòåìå ñîãëàñíî ïðîèçîøåäøåìó ñîáûòèþ
		reform(im, jm, vp1, vp2);

		retime(im);

		/*
		Åñëè jm > 0, òî ýòî çíà÷èò, ÷òî ïðîèçîøëî
		ñîóäàðåíèå ÷àñòèö â ñèñòåìå, óâåëè÷èâàåì
		ñ÷¸ò÷èê ñîóäàðåíèé
		*/
		if (jm >= 0) {
			++COLL_COUNT;
			retime(jm);
		}
	}

	/*
	Çàïóñêàåì ãëîáàëüíóþ ñèíõðîíèçàöèþ ÷àñòèö ïî ñîáñòâåííîìó âðåìåíè
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


void small_step(int devision) {
	particle p1;
	long int i, im, jm, vp1, vp2;
	double time = get_maximum_particle_time();

	COLL_COUNT = 0;
	jm = 0;

	while (COLL_COUNT < NP / devision || jm < 0) {
		// ñ÷èòûâàåì ïåðâîå ñîáûòèå èç ëèíåéêè ñîáûòèé
		im = time_queue[1].im;
		jm = time_queue[1].jm;
		vp1 = time_queue[1].vp1;
		vp2 = time_queue[1].vp2;

		// óäàëÿåì ïåðâîå ñîáûòèå
		delete_event(1);

		p1 = particles[im];

		/*
		Îáíóëÿåì óêàçàòåëü íà ñîáûòèå ÷àñòèö, ó÷àâñòâóþùèõ â ýòîì
		ñîáûòèè, ïðèðàâíèâàÿ èõ -1, òàê, ÷òîáû îíè íå óêàçûâàëè íà
		ñóùåñòâóþùèå ñîáûòèÿ.
		*/
		p1.ti = -1;
		if (jm >= 0) particles[jm].ti = -1;

		/*
		Åñëè íàñòóïèëî ñîáûòèå -100, çíà÷èò ïðîñòî
		ñèíõðîíèçèðóåì ÷àñòèöó ñ êîòîðîé ïðîèçîøëî
		äàííîå ñîáûòèå ñ ãëîáàëüíûì âðåìåíåì ñèñòåìû.
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

		dP_time += p1.dt;

		p1.dt = 0.0;
		time = p1.t;


		particles[im] = p1;

		// Ïðîèçâîäèì èçìåíåíèÿ â ñèñòåìå ñîãëàñíî ïðîèçîøåäøåìó ñîáûòèþ
		reform(im, jm, vp1, vp2);

		retime(im);

		/*
		Åñëè jm > 0, òî ýòî çíà÷èò, ÷òî ïðîèçîøëî
		ñîóäàðåíèå ÷àñòèö â ñèñòåìå, óâåëè÷èâàåì
		ñ÷¸ò÷èê ñîóäàðåíèé
		*/
		if (jm >= 0) {
			++COLL_COUNT;
			retime(jm);
		}
	}

	/*
	Çàïóñêàåì ãëîáàëüíóþ ñèíõðîíèçàöèþ ÷àñòèö ïî ñîáñòâåííîìó âðåìåíè
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
Ôóíêöèÿ ïîëó÷åíèÿ ïðîôèëÿ ïëîòíîñòè ñèñòåìû

Ïåðåä ñíÿòèåì õàðàêòåðèñòèê ïðîèçâîäèòñÿ ñèíõðîíèçàöèÿ âñåõ ÷àñòèö
ñ ãëîáàëüíûì âðåìåíåì ñèñòåìû.

Àðãóìåíòû:
steps - ÷èñëî çàìåðîâ ïëîòíîñòè äëÿ ïîñòðîåíèÿ äàííûõ î ïðîôèëå ïëîòíîñòè
accuracy - ÷èñëî ÿ÷ååê ïî Õ íà êîòîðûå ðàçáèâàåòñÿ ïðîôèëü ïëîòíîñòè
è â êîòîðûõ ñ÷èòàåòñÿ êîëè÷åñòâî ÷àñòèö
file_name - èìÿ ôàéëà äëÿ ñîõðàíåíèÿ äàííûõ
*/
void image(long steps, long accuracy, std::string file_name) {
	// ÷èñëî ÿ÷ååê, â êîòîðûõ áóäåì ðàññ÷èòûâàòü êîëè÷åñòâî ÷àñòèö
	const unsigned W = unsigned(L*accuracy);

	printf("\n W = %ld \n", W);

	// ìàññèâ, â êîòîðûé áóäåì ñîáèðàòü äàííûå ïî êîëè÷åñòâó ÷àñòèö
	long *img = new long[W];
	double t_global, dt;

	printf("INFO: Image started for %ld steps with accuracy %ld small cells per 1.0\n", steps, accuracy);

	// çàïîëíÿåì ìàññèâ íóëÿìè, ïðåæäå ÷åì ñ÷èòàòü êîëè÷åñòâî ÷àñòèö
	for (unsigned g = 0; g < W; g++) img[g] = 0;

	for (long h = 0; h < steps; h++) {
		step();  // äàëåì îäíî ñîóäàðåíèå íà ÷àñòèöó ìåæäó çàìåðàìè

				 // î÷èñòêà ýêðàíà è âûâîä èíôîðìàöèè î ïðîãðåññå
		printf("\r%76c\r", ' ');
		printf("%ld / %ld", long(h + 1), steps);

		/*
		Ðàññ÷èòûâàåì ãëîáàëüíîå âðåìÿ ñèñòåìû, êîòîðîå ðàâíî
		ñàìîìó áîëüøîìó ñîáñòâåííîìó âðåìåíè ÷àñòèö â ñèñòåìå.
		*/
		t_global = get_maximum_particle_time();

		for (long i = 0; i < NP; ++i) {
			/*
			Ñèíõðîíèçèðóåì ÷àñòèöó ñ ãëîáàëüíûì âðåìåíåì ñèñòåìû
			*/
			dt = t_global - particles[i].t;

			/*
			Îïðåäåëÿåì â êàêóþ èç ÿ÷ååê ïî X ïîïàäàåò äàííàÿ ÷àñòèöà
			ïîñëå ñèíõðîíèçàöèè å¸ ïî âðåìåíè ñ ãëîáàëüíûì âðåìåíåì ñèñòåìû
			*/
			long m = long((particles[i].x[0] + particles[i].vx * dt) * 10000.0 * accuracy) / 10000;
			img[m]++;  // óâåëè÷èâàåì êîëè÷åñòâî ÷àñòèö â ÿ÷åéêå íà 1.
		}
	}

	double x = 0.0, delta_x = 1.0 / accuracy;
	FILE *profile_file = fopen(file_name.c_str(), "w+");
	for (unsigned f = 0; f < W; f++) {
		/*
		Ñîõðàíÿåì íà÷àëüíóþ êîîðäèíàòó X äëÿ ÿ÷åéêè è ÷èñëî ÷àñòèö,
		êîòîðûå íàõîäèëèñü â äàííîì "ñëîå" âî âðåìÿ èçìåðåíèé
		*/
		fprintf(profile_file, "%f %ld\n", x, img[f]);
		x += delta_x;
	}
	fclose(profile_file);

	delete[] img;

	printf("\n INFO: Image completed. Information saved to file: %s \n", file_name.c_str());

	save("new.txt");
	load_seed("new.txt");
}


void function_g(long steps_count, double x1, double x2, std::string file_name) {

	long Ni = 0;
	double g_data[1000];
	double dt;

	for (int i = 0; i < 1000; i++)
		g_data[i] = 0;

	for (long w = 0; w < steps_count; w++) {
		step();

		//check_particles();

		// î÷èñòèòü ñòðîêó
		printf("\r%76c\r", ' ');
		printf("%ld / %ld", w + 1, steps_count);

		double t_global = get_maximum_particle_time();

		for (long i = 0; i < NP; i++) {
			particle p1 = particles[i];

			dt = t_global - p1.t;
			p1.x[0] = p1.x[0] + p1.vx * dt;
			p1.y[0] = p1.y[0] + p1.vy * dt;
			p1.z[0] = p1.z[0] + p1.vz * dt;

			if ((p1.x[0] >= x1) && (p1.x[0] <= x2) &&
				(p1.y[0] < A - 5.0) && (p1.y[0] > 5.0) &&
				(p1.z[0] < A - 5.0) && (p1.z[0] > 5.0))
			{
				Ni += 1;

				for (long j = 0; j < NP; j++) {
					if (i == j) continue;

					particle p2 = particles[j];

					dt = t_global - p2.t;
					p2.x[0] = p2.x[0] + p2.vx * dt;
					p2.y[0] = p2.y[0] + p2.vy * dt;
					p2.z[0] = p2.z[0] + p2.vz * dt;

					double dx = p1.x[0] - p2.x[0];
					double dy = p1.y[0] - p2.y[0];
					double dz = p1.z[0] - p2.z[0];
					double r = sqrt(dx*dx + dy*dy + dz*dz);

					if (r < 10.0L) {
						long u = long(r * 100);

						g_data[u] = g_data[u] + 1;
					}
				}
			}
		}
	}

	FILE *data_file = fopen(file_name.c_str(), "w+");
	fprintf(data_file, "%d\n", steps_count);

	double r = 0.0;
	for (int i = 0; i < 1000; i++) {
		g_data[i] = g_data[i] / (4.0 * PI * ((r + 0.01)*(r + 0.01)*(r + 0.01) - r*r*r) / 3.0);
		fprintf(data_file, "%.5le\n", g_data[i] / Ni);

		r += 0.01;
	}
	fclose(data_file);

}


/*
void count_of_nearest_particles(int steps_count, double x1, double x2, std::string file_name) {

double Ni = 0, count_of_particles = 0;

for (int w = 0; w < steps_count; w++) {
step();

// î÷èñòèòü ñòðîêó
printf("\r%76c\r", ' ');
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



void profile_medium(double x1, double x2, int dots_per_particle, long steps, std::string file_name) {
	double x, dt;
	double *y = new double[1000];
	double *z = new double[1000];
	long *f = new long[1000];
	long i, j, w = 0, u, m, particles_count;

	for (i = 0; i < 1000; i++) {
		y[i] = 0.0;
		z[i] = 0.0;
		f[i] = 0;
	}

	// îòêðûâàåì ôàéë íà çàïèñü
	FILE *profile_file = fopen(file_name.c_str(), "w+");

	printf("INFO: Profile started: %.5le - %.5le.\n", x1, x2);

	double t_global = get_maximum_particle_time();

	particles_count = 0;
	for (j = 0; j < NP; j++) {
		/*
		Ñèíõðîíèçèðóåì ÷àñòèöó ïî âðåìåíè ñ ãëîáàëüíûì
		âðåìåíåì ñèñòåìû è ïðîâåðÿåì å¸ êîîðäèíàòó Õ,
		åñëè ÷àñòèöà íàõîäèòñÿ â óêàçàííîì äèàïàçîíå,
		òî ñîõðàíÿåì å¸ Y è Z êîîðäèíàòû â ôàéë.
		*/
		dt = t_global - particles[j].t;
		x = particles[j].x[0] + particles[j].vx * dt;

		if (x <= x2 && x >= x1) {
			f[particles_count] = j;
			particles_count++;
		}
	}

	printf("\n Particles count in this layer: %ld\n", particles_count);

	fprintf(profile_file, "%ld\n", particles_count);

	for (i = 0; i < dots_per_particle; i++) {

		/*
		Ìåæäó èçìåðåíèÿìè äàííûõ î ïîëîæåíèÿõ ÷àñòèö
		äåëàåì óêàçàííîå êîëè÷åñòâî ñîóäàðåíèé íà ÷àñòèöó â ñèñòåìå
		*/
		for (u = 0; u < steps; u++) {
			step();

			// î÷èñòèòü ñòðîêó
			printf("\r%76c\r", ' ');
			printf("%ld / %ld", i*steps + u, steps*dots_per_particle);

			t_global = get_maximum_particle_time();

			w = 0;
			for (j = 0; j < particles_count; j++) {
				/*
				Ñèíõðîíèçèðóåì ÷àñòèöó ïî âðåìåíè ñ ãëîáàëüíûì
				âðåìåíåì ñèñòåìû è ïðîâåðÿåì å¸ êîîðäèíàòó Õ,
				åñëè ÷àñòèöà íàõîäèòñÿ â óêàçàííîì äèàïàçîíå,
				òî ñîõðàíÿåì å¸ Y è Z êîîðäèíàòû â ôàéë.
				*/
				dt = t_global - particles[f[j]].t;
				x = particles[f[j]].x[0] + particles[f[j]].vx * dt;

				y[w] += particles[f[j]].y[0] + particles[f[j]].vy * dt;
				z[w] += particles[f[j]].z[0] + particles[f[j]].vz * dt;

				w++;
			}
		}

		for (m = 0; m < w; m++) {
			// ñîõðàíÿåì Y è Z êîîðäèíàòû ÷àñòèöû
			fprintf(profile_file, "%ld\n", f[m]);
			fprintf(profile_file, "%.15le\n", y[m] / steps);
			fprintf(profile_file, "%.15le\n", z[m] / steps);

			y[m] = 0;
			z[m] = 0;
		}
	}

	fclose(profile_file);

	printf("INFO: Profile completed. Information saved to file: %s\n", file_name.c_str());

	delete[] y;
	delete[] z;
	delete[] f;

	save("new.txt");
	load_seed("new.txt");
}

/*
Ôóíêöèÿ ñîçäàíèÿ ðàçðåçà ñèñòåìû â ïèêå,
ïîçâîëÿåò ïðîñìàòðèâàòü ðàñïîëîæåíèå ÷àñòèö â ñëîå,
ïàðàëëåëüíîì èäåàëüíîé ñòåíêå.

Àðãóìåíòû:
x1 - íà÷àëüíàÿ X êîîðäèíàòà "ñðåçà"
x2 - êîíå÷íàÿ X êîîðäèíàòà "ñðåçà"
dots_per_particle - êîëè÷åñòâî ñíÿòèé äàííûõ ÷åðåç ðàâíîå êîëè÷åñòâî ñîóäàðåíèé â ñèñòåìå
steps - ÷èñëî ñîóäàðåíèé ìåæäó äâóìÿ ñíÿòèÿìè äàííûõ î ïîëîæåíèè ÷àñòèö â âûáðàííîì ñëîå
file_name - èìÿ ôàéëà äëÿ ñîõðàíåíèÿ äàííûõ
*/
void profile(double x1, double x2, long int dots_per_particle, long int steps, std::string file_name) {
	double x, y, z, dt;
	long int i, j;

	// îòêðûâàåì ôàéë íà çàïèñü
	FILE *profile_file = fopen(file_name.c_str(), "w+");

	printf("INFO: Profile started: %.5le - %.5le.\n", x1, x2);

	for (i = 0; i < dots_per_particle; ++i) {

		// î÷èñòèòü ñòðîêó
		printf("\r%76c\r", ' ');
		printf("%ld / %ld", i, dots_per_particle);

		/*
		Ìåæäó èçìåðåíèÿìè äàííûõ î ïîëîæåíèÿõ ÷àñòèö
		äåëàåì óêàçàííîå êîëè÷åñòâî ñîóäàðåíèé íà ÷àñòèöó â ñèñòåìå
		*/
		for (j = 0; j < steps; j++)
			step();

		double t_global = get_maximum_particle_time();

		for (j = 0; j < NP; j++) {
			/*
			Ñèíõðîíèçèðóåì ÷àñòèöó ïî âðåìåíè ñ ãëîáàëüíûì
			âðåìåíåì ñèñòåìû è ïðîâåðÿåì å¸ êîîðäèíàòó Õ,
			åñëè ÷àñòèöà íàõîäèòñÿ â óêàçàííîì äèàïàçîíå,
			òî ñîõðàíÿåì å¸ Y è Z êîîðäèíàòû â ôàéë.
			*/
			dt = t_global - particles[j].t;
			x = particles[j].x[0] + particles[j].vx * dt;

			if (x <= x2 && x >= x1) {
				y = particles[j].y[0] + particles[j].vy * dt;
				z = particles[j].z[0] + particles[j].vz * dt;

				// ñîõðàíÿåì Y è Z êîîðäèíàòû ÷àñòèöû
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
x1, x2 - X êîîðäèíàòû äâóõ ïëîñêîñòåé, îãðàíè÷èâàþùèõ ñëîé êðèñòàëëà
m1 - êîëè÷åñòâî èçìåðåíèé, êîòîðûå áóäó ñäåëàíû â òå÷åíèè îäíîé èòåðàöèè
m2 - êîëè÷åñòâî èçìåðåíèé, êîòîðûå äåëàþòñÿ äëÿ íàõîæäåíèÿ ñðåäíåãî ïîëîæåíèÿ öåíòðà ÷àñòèöû
m3 - êîëè÷åñòâî èòåðàöèé èçìåðåíèé
m4 - êîëè÷åñòâî ñîóäàðåíèé â ñðåäíåì íà îäíó ÷àñòèöó â ñèñòåìå, êîòîðûå äîëæíû ïðîèîéòè
ìåæäó äâóìÿ áëèæàéøèìè èçìåðåíèÿìè.
*/
void F1(double x1, double x2, int m1, int m2, int m3, int m4, std::string file_name) {

	double t_global, dt;
	particle p1;

	FILE *f1_file = fopen(file_name.c_str(), "w+");
	FILE *debug_file = fopen("debug.log", "w+");

	// Äèíàìè÷åñêè ñîçäà¸ì äâóìåðíûå ìàññèâû (ýòî íåîáõîäèìî äåëàòü äèíàìè÷åñêè,
	// òàê êàê îíè áîëüøîãî ðàçìåðà è ìû íå ìîæåì ñîçäàòü èõ â ñòåêå,
	// ê òîìó æå, ðàçìåð ìàññèâà çàâèñèò îò çíà÷åíèé ïåðåìåííîé m1)
	// 200 - ïî êîëè÷åñòâó ÷àñòèö â ñëîå, ñ çàïàñîì (÷àñòèö áóäåò 100+, íî ìåíüøå 200)
	double **x = new double*[200];
	double **y = new double*[200];
	double **z = new double*[200];
	for (int i = 0; i < 200; i++) {
		x[i] = new double[m1 + 2];
		y[i] = new double[m1 + 2];
		z[i] = new double[m1 + 2];
	}

	double *middle_x1 = new double[200];
	double *middle_y1 = new double[200];
	double *middle_z1 = new double[200];
	double *middle_x2 = new double[200];
	double *middle_y2 = new double[200];
	double *middle_z2 = new double[200];

	long *list_of_particles = new long[200];
	long *list_of_particles_alfa = new long[50];
	long *list_of_particles_betta = new long[50];

	long particles_count, particles_count_alfa, particles_count_betta, w, w1, w2;

	for (int iter = 0; iter < m3; iter++) {

		// ×èñëî ÷àñòèö, êîòîðûå áóäóò âûáðàíû äëÿ íàáëþäåíèÿ â äàííîì ñëîå:
		particles_count = 0;

		// ×èñëî ÷àñòèö, öåíòðû êîòîðûõ ðàñïîëîæåíû âäîëü îäíîé èç îñåé ñèììåòðèè ñëîÿ:
		particles_count_alfa = 0;
		// ×èñëî ÷àñòèö, öåíòðû êîòîðûõ ðàñïîëîæåíû âäîëü åùå îäíîé îñè ñèììåòðèè ñëîÿ
		// äàííàÿ îñü ïåðåñåêàåòñÿ ñ ïåðâîé â öåíòðå îäíîé èç ÷àñòèö):
		particles_count_betta = 0;

		// Çàïîëíÿåì âñå ìàññèâû íóëÿìè:
		for (w1 = 0; w1 < 200; w1++) {
			middle_x1[w1] = 0.0;
			middle_y1[w1] = 0.0;
			middle_z1[w1] = 0.0;
			middle_x2[w1] = 0.0;
			middle_y2[w1] = 0.0;
			middle_z2[w1] = 0.0;

			list_of_particles[w1] = -1;

			for (w2 = 0; w2 < m1 + 1; w2++) {
				x[w1][w2] = 0.0;
				y[w1][w2] = 0.0;
				z[w1][w2] = 0.0;
			}
		}

		// Íå ðàññìàòðèâàåì ÷àñòèöû, íàõîäÿùèåñÿ âáëèçè ïåðèîäè÷åñêèõ
		// ãðàíè÷íûõ óñëîâèé, à òàê æå ðàññìàòðèâàåì òîëüêî ÷àñòèöû,
		// íàõîäÿùèåñÿ â äàííîì ñëîå.
		// Â ýòîì öèêëå ìû ñíà÷àëà ïåðåìåùàåì âñå ÷àñòèöû âî âðåìÿ ñèñòåìû,
		// à çàòåì ïðîâåðÿåì óäîâëåòâîðÿåò ëè òåêóùåå ïîëîæåíèå öåíòðà ÷àñòèöû êðèòåðèÿì,
		// ÷òîáû íàéòè âñå ÷àñòèöû, ïîëîæåíèå êîòîðûõ íåîáõîäèìî îòñëåæèâàòü:
		t_global = get_maximum_particle_time();
		for (long i = 0; i < NP; i++) {
			p1 = particles[i];
			dt = t_global - p1.t;
			p1.x[0] = p1.x[0] + p1.vx * dt;
			p1.y[0] = p1.y[0] + p1.vy * dt;
			p1.z[0] = p1.z[0] + p1.vz * dt;

			// Îòáèðàåì òîëüêî ÷àñòèöû â äàííîì ñëîå òàêèå, êîòîðûå íàõîäÿòñÿ íà ðàññòîÿíèè
			// áîëüøå 6 ðàäèóñîâ îò ïåðèîäè÷åñêèõ ãðàíè÷íûõ óñëîâèé:
			if (p1.x[0] >= x1 && p1.x[0] <= x2 &&
				p1.y[0] >= 6.0 * particle_R && p1.y[0] <= A - 6.0 * particle_R &&
				p1.z[0] >= 6.0 * particle_R && p1.z[0] <= A - 6.0 * particle_R) {
				list_of_particles[particles_count] = i;
				particles_count++;
			}
		}

		for (long i = 0; i < m1; i++) {
			// Âûâîäèì íà ýêðàí òåêóùèé ïðîãðåññ ðàñ÷åòà:
			printf("\r%76c\r", ' ');
			printf("%ld / %ld", iter*m1 + i + 1, m1*m3);

			// Äåëàåì íåñêîëüêî ñîóäàðåíèé íà ÷àñòèöó ìåæäó èçìåðåíèÿìè
			// ÷òîáû ïîëîæåíèÿ ÷àñòèö áûëè äîñòàòî÷íî ñëó÷àéíû.
			for (int r = 0; r < m4; r++) {
				step();
			}

			t_global = get_maximum_particle_time();
			for (int m = 0; m < particles_count; m++) {
				w = list_of_particles[m];

				// Cèíõðîíèçèðóåì âðåìÿ ÷àñòèöû ñ âðåìåíåì ñèñòåìû
				p1 = particles[w];
				dt = t_global - p1.t;
				p1.x[0] = p1.x[0] + p1.vx * dt;
				p1.y[0] = p1.y[0] + p1.vy * dt;
				p1.z[0] = p1.z[0] + p1.vz * dt;

				// Cóììèðóåì çíà÷åíèÿ êîîðäèíàò öåíòðîâ óçëîâ çà ïåðâûå m2 è çà 
				// ïîñëåäíèå m2 ñîóäàðåíèé
				if (i < m2) {
					middle_x1[m] += p1.x[0];
					middle_y1[m] += p1.y[0];
					middle_z1[m] += p1.z[0];
				}
				else if (i >= m1 - m2) {
					middle_x2[m] += p1.x[0];
					middle_y2[m] += p1.y[0];
					middle_z2[m] += p1.z[0];
				}

				// Çàïîìèíàåì âñå êîîðäèíàòû âñåõ ÷àñòèö
				x[m][i] = p1.x[0];
				y[m][i] = p1.y[0];
				z[m][i] = p1.z[0];
			}
		}

		// Íàõîäèì ñðåäíèå çíà÷åíèÿ êîîðäèíàò öåíòðîâ óçëîâ çà ïåðâûå m2 è çà 
		// ïîñëåäíèå m2 ñîóäàðåíèé
		for (int m = 0; m < particles_count; m++) {
			middle_x1[m] /= double(m2);
			middle_y1[m] /= double(m2);
			middle_z1[m] /= double(m2);
			middle_x2[m] /= double(m2);
			middle_y2[m] /= double(m2);
			middle_z2[m] /= double(m2);
		}

		long min_m, repeat;
		double min_y = A, min_z = A, dy, dz, alfa, alfa2 = -1.0, betta, betta2 = -1.0;

		// Íà÷èíàåì ïîèñê ëèíèè, ïî êîòîðîé áóäåì îïðåäåëÿòü ïîâîðîò 
		// îñåé êîîðäèíàò äëÿ èññëåäóåìîãî ñëîÿ

		// Äëÿ íà÷àëà íàõîäèì ñàìóþ ëåâóþ íèæíþþ ÷àñòèöó â ñëîå, èç ÷èñëà òåõ,
		// êîòîðûå áûëè îòîáðàíû äëÿ ðàññìîòðåíèÿ (íå áëèæå ÷åì 6 ðàäèóñîâ îò
		// ïåðèîäè÷åñêèõ ãðàíè÷íûõ óñëîâèé):
		for (int m = 0; m < particles_count; m++)
		{
			if (((middle_y1[m] + middle_y2[m]) / 2.0 < min_y) &&
				((middle_z1[m] + middle_z2[m]) / 2.0 < min_z))
			{
				min_y = middle_y1[m];
				min_z = middle_z1[m];
				min_m = m;
			}
		}

		list_of_particles_alfa[0] = min_m;
		particles_count_alfa = 1;
		repeat = 1;
		while (repeat > 0)
		{
			repeat = -1;
			for (int m = 0; repeat < 0 && m < particles_count; m++) {
				if (m != min_m) {
					dy = middle_y1[m] - middle_y1[min_m];
					dz = middle_z1[m] - middle_z1[min_m];

					// Èùåì áëèæàéøóþ ê äàííîé ÷àñòèöå ÷àñòèöó íà ðàññòîÿíèè íå äàëüøå 3.2 ðàäèóñîâ
					// òàê, ÷òîáû ëèíèÿ, ïåðåñåêàþùàÿ öåíòðû èõ óçëîâ,
					// îáðàçîâûâàëà îñòðûé óãîë ñ îñüþ OY.
					if ((sqrt(dy*dy + dz*dz) < 3.2*particle_R) && (dy > 0.0 && dz > 0.0))
					{
						alfa = atan2(dy, dz);

						// Åñëè òàêàÿ ÷àñòèöà íàéäåíà è îíà ëåæèò íà òîé æå ïðÿìîé, ÷òî
						// è âñå óæå îòîáðàííûå äëÿ ýòîé ëèíèè ðàíåå ÷àñòèöû, òî çàïèñûâàåì åå â ñïèñîê:
						if (alfa2 < 0.0 || fabs(alfa - alfa2) < 0.5)
						{
							list_of_particles_alfa[particles_count_alfa] = m;
							particles_count_alfa++;
							min_m = m;
							alfa2 = alfa;
							repeat = 1;
						}
					}
				}
			}
		}

		min_m = list_of_particles_alfa[particles_count_alfa / 2];
		list_of_particles_betta[0] = min_m;
		particles_count_betta = 1;
		repeat = 1;
		while (repeat > 0)
		{
			repeat = -1;
			for (int m = 0; repeat < 0 && m < particles_count; m++) {
				if (m != min_m) {
					dy = middle_y1[m] - middle_y1[min_m];
					dz = middle_z1[m] - middle_z1[min_m];

					// Èùåì ÷àñòèöó íà ðàññòîÿíèè íå äàëüøå 3.2 ðàäèóñîâ
					// òàê, ÷òîáû ëèíèÿ, ïåðåñåêàþùàÿ öåíòðû óçëîâ,
					// îáðàçîâûâàëà îñòðûé óãîë ñ îñüþ OY.
					if ((sqrt(dy*dy + dz*dz) < 3.2*particle_R) && (dy > 0.0 && dz < 0.0))
					{
						betta = atan2(dy, dz);

						if (betta2 < 0.0 || fabs(betta - betta2) < 0.5) {
							list_of_particles_betta[particles_count_betta] = m;
							particles_count_betta++;
							min_m = m;
							betta2 = betta;
							repeat = 1;
						}
					}
				}
			}
		}

		// Áåðåì ñðåäíþþ ÷àñòèöó â ñïèñêå ÷àñòèö âäîëü îäíîé èç îñåé ñèììåòðèè ñëîÿ
		// è íàõîäèì äëÿ ýòîé ÷àñòèöû äðóãóþ îñü ñèììåòðèè, ñîñòàâëÿÿ ñïèñîê âñåõ ÷àñòèö,
		// öåíòðû êîòîðûõ íàõîäÿòñÿ âäîëü ýòîé ëèíèè:
		min_m = list_of_particles_alfa[particles_count_alfa / 2];
		repeat = 1;
		while (repeat > 0)
		{
			repeat = -1;
			for (int m = 0; repeat < 0 && m < particles_count; m++) {
				if (m != min_m) {
					dy = middle_y1[min_m] - middle_y1[m];
					dz = middle_z1[min_m] - middle_z1[m];

					// Èùåì ÷àñòèöó íà ðàññòîÿíèè íå äàëüøå 3.2 ðàäèóñîâ
					// òàê, ÷òîáû ëèíèÿ, ïåðåñåêàþùàÿ öåíòðû óçëîâ,
					// îáðàçîâûâàëà îñòðûé óãîë ñ îñüþ OY.
					if ((sqrt(dy*dy + dz*dz) < 3.2*particle_R) && (dy > 0.0 && dz < 0.0))
					{
						betta = atan2(dy, dz);

						if (betta2 < 0.0 || fabs(betta - betta2) < 0.5) {
							list_of_particles_betta[particles_count_betta] = m;
							particles_count_betta++;
							min_m = m;
							betta2 = betta;
							repeat = 1;
						}
					}
				}
			}
		}

		double B_alfa1 = 0.0, B_alfa2 = 0.0, B_betta1 = 0.0, B_betta2 = 0.0;
		double K_alfa1 = 0.0, K_alfa2 = 0.0, K_betta1 = 0.0, K_betta2 = 0.0;
		double sum_y = 0.0, sum_z = 0.0, sum_yz = 0.0, sum_zz = 0.0;
		double sum_y2 = 0.0, sum_z2 = 0.0, sum_yz2 = 0.0, sum_zz2 = 0.0;

		// Ïîëó÷àåì óãîë ïîâîðîòà ëèíèè ìåòîäîì íàèìåíüøèõ êâàäðàòîâ
		// äëÿ âñåõ ÷àñòèö â ìàññèâå list_of_particles_alfa (ýòî îäíà èç
		// îñåé êðèñòàëëà â äàííîì ñëîå, ýòà îñü îáðàçóåò îñòðûé óãîë ñ îñüþ OY
		// è ïðîõîäèò ÷åðåç 10+ óçëîâ â äàííîì ñëîå):
		for (int m = 0; m < particles_count_alfa; m++) {
			w = list_of_particles_alfa[m];

			sum_y += middle_y1[w];
			sum_z += middle_z1[w];
			sum_yz += middle_y1[w] * middle_z1[w];
			sum_zz += middle_z1[w] * middle_z1[w];

			sum_y2 += middle_y2[w];
			sum_z2 += middle_z2[w];
			sum_yz2 += middle_y2[w] * middle_z2[w];
			sum_zz2 += middle_z2[w] * middle_z2[w];
		}

		// K_alfa1 - ýòî êîýôôèöèåíò â óðàâíåíèè ïðÿìîé (z = ky + b), ðàññ÷èòàííûé íà îñíîâàíèè
		// ïîëîæåíèÿ ÷àñòèö â äàííîì ñëîå íà ïðîòÿæåíèè ïåðâûõ m2 ñîóäàðåíèé:
		K_alfa1 = (particles_count_alfa*sum_yz - sum_y*sum_z) / (particles_count_alfa*sum_zz - sum_z*sum_z);
		// K_alfa1 - ýòî êîýôôèöèåíò â óðàâíåíèè ïðÿìîé (z = ky + b), ðàññ÷èòàííûé íà îñíîâàíèè
		// ïîëîæåíèÿ ÷àñòèö â äàííîì ñëîå íà ïðîòÿæåíèè ïîñëåäíèõ m2 ñîóäàðåíèé:
		K_alfa2 = (particles_count_alfa*sum_yz2 - sum_y2*sum_z2) / (particles_count_alfa*sum_zz2 - sum_z2*sum_z2);

		alfa = atan(K_alfa1);
		alfa2 = atan(K_alfa2);

		B_alfa1 = (sum_y - K_alfa1*sum_z) / particles_count_alfa;
		B_alfa2 = (sum_y - K_alfa2*sum_z) / particles_count_alfa;

		sum_y = 0.0; sum_z = 0.0; sum_yz = 0.0; sum_zz = 0.0;
		sum_y2 = 0.0; sum_z2 = 0.0; sum_yz2 = 0.0; sum_zz2 = 0.0;

		// Ïîëó÷àåì óãîë ïîâîðîòà âòîðîé ëèíèè ìåòîäîì íàèìåíüøèõ êâàäðàòîâ
		// (äàííàÿ ëèíèÿ íå ñîâïàäàåò ñ ïåðâîé ëèíèåé, íî îíà ïåðåñåêàåòñÿ ñ íåé
		// â òî÷êå, ñîâïàäàþùåé ñ öåíòðîì îäíîãî èç óçëîâ, ÷åðåç êîòîðûå ïðîõîäèò
		// ïåðâàÿ ëèíèÿ):
		for (int m = 0; m < particles_count_betta; m++) {
			w = list_of_particles_betta[m];

			sum_y += middle_y1[w];
			sum_z += middle_z1[w];
			sum_yz += middle_y1[w] * middle_z1[w];
			sum_zz += middle_z1[w] * middle_z1[w];

			sum_y2 += middle_y2[w];
			sum_z2 += middle_z2[w];
			sum_yz2 += middle_y2[w] * middle_z2[w];
			sum_zz2 += middle_z2[w] * middle_z2[w];
		}

		K_betta1 = (particles_count_betta*sum_yz - sum_y*sum_z) / (particles_count_betta*sum_zz - sum_z*sum_z);
		K_betta2 = (particles_count_betta*sum_yz2 - sum_y2*sum_z2) / (particles_count_betta*sum_zz2 - sum_z2*sum_z2);

		B_betta1 = (sum_y - K_betta1*sum_z) / particles_count_betta;
		B_betta2 = (sum_y - K_betta2*sum_z) / particles_count_betta;

		double delta_alfa = (alfa2 - alfa) / (m1 - m2); // Íàõîäèì ñêîðîñòü âðàùåíèÿ ñëîÿ (â ãðàäóñàõ çà âðåìÿ â ñðåäíåì ìåæäó
														// äâóìÿ èçìåðåíèÿìè)
		fprintf(debug_file, "\n /////////////////////////////////////////////////////\n");
		fprintf(debug_file, "\n delta alfa = %.5le\n", delta_alfa);

		// Íàõîäèì ïåðåñå÷åíèå äâóõ íàéäåííûõ îñåé êðèñòàëëà â äàííîì ñëîå
		// ÷òîáû ýòî ñäåëàòü - ðåøàåì óðàâíåíèå ïåðåñå÷åíèÿ äâóõ ïðÿìûõ, îáðàçîâàííûõ
		// âûáðàííûìè íà ïðåäûäóùèõ øàãàõ ÷àñòèö:
		// y1, z1 - êîîðäèíàòû òî÷êè ïåðåñå÷åíèÿ ýòèõ ëèíèé çà ïåðâûå 200 ñîóäàðåíèé,
		// y2, z2 - êîîðäèíàòû òî÷êè ïåðåñå÷åíèÿ ëèíèé çà ïîñëåäíèå 200 ñîóäàðåíèé.
		double y1 = K_alfa1*(B_betta1 - B_alfa1) / (K_alfa1 - K_betta1) + B_alfa1;
		double y2 = K_alfa2*(B_betta2 - B_alfa2) / (K_alfa2 - K_betta2) + B_alfa2;
		double z1 = (B_betta1 - B_alfa1) / (K_alfa1 - K_betta1);
		double z2 = (B_betta2 - B_alfa2) / (K_alfa2 - K_betta2);

		fprintf(debug_file, "\n delta Y = %.5le \n delta Z = %.5le", y2 - y1, z2 - z1);

		fprintf(debug_file, "\n Y1 = %.5le, Z1 = %.5le \n", y1, z1);
		min_m = list_of_particles_alfa[particles_count_alfa / 2];
		fprintf(debug_file, "\n min_m Y = %.5le, Z = %.5le \n", middle_y1[min_m], middle_z1[min_m]);

		// Ðàññ÷èòûâàåì ñìåùåíèå âñåãî ñëîÿ êðèñòàëëà:
		double delta_y = (y2 - y1) / (m1 - m2);
		double delta_z = (z2 - z1) / (m1 - m2);

		int k1 = 1;
		int k2 = 1;
		for (int j = 0; j < particles_count; j++) {

			if (abs(middle_y2[j] - middle_y1[j]) > abs(y1 - y2)) {
				w = list_of_particles[j];
				fprintf(debug_file, "\n %d Particle %d, dY = %.5le", k1, w, middle_y2[j] - middle_y1[j]);
				k1++;
			}

			if (abs(middle_z2[j] - middle_z1[j]) > abs(z1 - z2)) {
				w = list_of_particles[j];
				fprintf(debug_file, "\n %d Particle %d, dZ = %.5le", k2, w, middle_z2[j] - middle_z1[j]);
				k2++;
			}
		}

		FILE *debug_file2 = fopen("debug2.log", "w+");
		for (int j = 0; j < particles_count; j++) {
			fprintf(debug_file2, "%d %.5le %.5le \n", j, abs(middle_y1[j] - middle_y2[j]), abs(middle_z1[j] - middle_z2[j]));

			double y_min = A;
			double z_min = A;
			double y_max = 0;
			double z_max = 0;

			for (int i = m2 / 2; i < m1 - m2 / 2; i++) {
				if (y[j][i] < y_min) y_min = y[j][i];
				if (z[j][i] < z_min) z_min = z[j][i];
				if (y[j][i] > y_max) y_max = y[j][i];
				if (z[j][i] > z_max) z_max = z[j][i];
			}

			fprintf(debug_file2, "%d %.5le %.5le \n", j, abs(y_max - y_min), abs(z_max - z_min));
		}
		fclose(debug_file2);

		// Ïðîõîäèì ïî âñåì çàïîìíåííûì íàìè êîîðäèíàòàì ñ m2/2 ñîóäàðåíèÿ
		// íà ÷àñòèöó äî m1-m2/2 ñîóäàðåíèÿ, òî åñòü áåðåì äëÿ ðàññìîòðåíèÿ m1-m2 ñîóäàðåíèé:
		for (int i = m2 / 2; i < m1 - m2 / 2; i++) {
			// Äëÿ êàæäîé ÷àñòèöû â äàííîì ñëîå (åñëè ÷àñòèöà ïðîøëà âñå óñëîâèÿ îòáîðà):
			for (int j = 0; j < particles_count; j++) {

				// Çàïðàøèâàåì íîìåð ÷àñòèö â ñïèñêå ÷àñòèö äàííîãî ñëîÿ:
				w = list_of_particles[j];

				// Ðàññ÷èòûâàåì íîìåð òåêóùåé èòåðàöèè ÷òîáû ñäåëàòü ïîïðàâêó íà ïîâîðîò è
				// ñìåùåíèå ñëîÿ êðèñòàëëà:
				int f = i - m2 / 2;

				// Ðàññ÷èòûâàåì êîîðäèíàòû ÷àñòèöû ñ ó÷åòîì ñìåùåíèÿ ñëîÿ è ñðåäíåãî ïîëîæåíèÿ
				// öåíòðà äàííîé ÷àñòèöû (òî åñòü ðàññ÷èòûâàåì òåêóùåå ñìåùåíèå ÷àñòèöû íà äàííîì øàãå îò
				// ñðåäíåãî ïîëîæåíèÿ öåíòðà óçëà, êîòîðûé îáðàçîâàí ýòîé ÷àñòèöåé):
				double x_mol = x[j][i] - (middle_x1[j] + middle_x2[j]) / 2.0;
				double y_mol0 = (middle_y1[j] - y[j][i] + delta_y*f)*cos(alfa) + (middle_z1[j] - z[j][i] + delta_z*f)*sin(alfa);
				double z_mol0 = -(middle_y1[j] - y[j][i] + delta_y*f)*sin(alfa) + (middle_z1[j] - z[j][i] + delta_z*f)*cos(alfa);

				double y_mol = y_mol0*cos(delta_alfa*f) + z_mol0*sin(delta_alfa*f);
				double z_mol = -y_mol0*sin(delta_alfa*f) + z_mol0*cos(delta_alfa*f);

				// Çàïèñûâàåì íîìåð ÷àñòèöû è ðàññ÷èòàííîå ñìåùåíèå â ôàéë:
				fprintf(f1_file, "%d %.5le %.5le %.5le\n", w, x_mol, y_mol, z_mol);
			}
		}
	}

	// Çàêðûâàåì âñå ôàéëû, êîòîðûå ìû ðåäàêòèðîâàëè:
	fclose(f1_file);
	fclose(debug_file);

	// Óäàëÿåì èç ïàìÿòè âñå äèíàìè÷åñêèå ìàññèâû:
	for (int i = 0; i < 200; i++) {
		delete[] x[i];
		delete[] y[i];
		delete[] z[i];
	}
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] middle_x1;
	delete[] middle_y1;
	delete[] middle_z1;
	delete[] middle_x2;
	delete[] middle_y2;
	delete[] middle_z2;
	delete[] list_of_particles;
	delete[] list_of_particles_alfa;
}


/*
Ôóíêöèÿ ñæàòèÿ ñèñòåìû, ïîçâîëÿåò ñæèìàòü èëè ðàñøèðÿòü ñèñòåìó äî çàäàííîé ïëîòíîñòè
ñ çàäàííûì ìàêñèìàëüíûì øàãîì ïî ïëîòíîñòè.

Àðãóìåíòû:
compress_to_etta - èòîãîâàÿ ïëîòíîñòü
delta_max - ìàêñèìàëüíûé ðàçðåø¸ííûé øàã ïî L
KK1 - êîëè÷åñòâî ñîóäàðåíèé íà îäíó ÷àñòèöó â ñèñòåìå ìåæäó
ìàëåíüêèìè ñæàòèÿìè ïî ïëîòíîñòè
KK2 - êîëè÷åñòâî ñäâèãîâ ñòåíêè, ïîñëå êîòîðûõ íåîáõîäèìî ïðîâåñòè
äîïîëíèòåëüíî KK3 ñîóäàðåíèé
KK3 - êîëè÷åñòâî ñîóäàðåíèå íà îäíó ÷àñòèöó â ñèñòåìå êàæäûå KK2
ñäâèãîâ ñòåíêè
type - òèï ñæàòèÿ:
0 - ïîäîäâèãàòü òîëüêî ëåâóþ ñòåíêó
1 - ïîäîäâèãàòü òîëüêî ïðàâóþ ñòåíêó
2 - ïîäîäâèãàòü îäíîâðåìåííî îáå ñòåíêè íà îäèíàêîâîå ðàññòîÿíèå
*/

void compress(double compress_to_etta, double delta_max, long int KK1, long int KK2, long int KK3, int type) {
	long int compression_steps_done = 0;
	double etta = (4.0 * PI * particle_R3 * NP) / (3.0 * A * A * (L - 2.0*particle_R));
	double min = 1.0e+100, max = -1.0, x, deltaL, dx;
	double t_global, dt;
	double L_ideal = ((4.0 * PI * particle_R3 * NP) / compress_to_etta) / (3.0 * A * A) + 2.0*particle_R;
	double p, q, Qf, R_ideal = 0.0;

	printf("\n INFO: Start to change system density... \n");

	if (type < 3) {
		printf(" \n Maximum delta_L = %.15le \n", delta_max);
	}
	else {
		printf(" \n Maximum delta_R = %.15le \n", delta_max);
		// Íàõîäèì èäåàëüíûé ðàäèóñ äëÿ ÷àñòèö, ðåøàÿ óðàâíåíèå
		// ìåòîäîì Êàðäàíî:
		p = (1.5 * compress_to_etta * A * A) / (PI * NP);
		q = -(0.75 * compress_to_etta * A * A * L) / (PI * NP);
		Qf = sqrt(q*q / 4.0 + p*p*p / 27.0);
		R_ideal = pow(-q / 2.0 + Qf, 1.0 / 3.0) + pow(-q / 2.0 + Qf, 1.0 / 3.0);
		printf("\n R_ideal %.15le \n", R_ideal);
	}

	printf(" Program will wait %ld collissions per particle between each changes in density.\n", KK1);
	printf(" Program will wait %ld collissions per particle between each %ld changes in density.\n", KK3, KK2);
	printf(" Type of compression: ");
	if (type == 0) printf("only position of left wall will be changed");
	if (type == 1) printf("only position of right wall will be changed");
	if (type == 2) printf("position of both walls will be changed");
	if (type == 3) printf("R of particles will be changed");
	printf("\n\n");

	// çàäà¸ì ïëîòíîñòü ñ òî÷íîñòüþ â 12 çíàêîâ
	while (abs(etta - compress_to_etta) > 1.0e-12) {
		if (type < 3) {
			if (etta < compress_to_etta) {
				max = -100.0;
				min = 1.1e+10;

				// Ñìåùàåì âñå ÷àñòèöû â òåêóùåå âðåìÿ ñèñòåìû ÷òîáû ñèíõðîíèçîâàòü
				// âñå ÷àñòèöû è çàòåì íàõîäèì êîîðäèíàòû ÷àñòèö, íàèáîëåå áëèçêî
				// ðàñïîëîæåííûõ îò èäåàëüíûõ ñòåíîê, ÷òîáû çíàòü íà ñêîëüêî ìîæíî
				// ïîäîäâèíóòü ñòåíêó íå êîñíóâøèñü ÷àñòèö.
				// Òàê æå îïðåäåëÿåì ìèíèìàëüíîå ðàñòîÿíèå íà êîòîðîì ÷àñòèöû íàõîäÿòñÿ
				// äðóã îò äðóãà, åñëè èçìåíåíèå ïëîòíîñòè ñèñòåìû ïðîèçâîäèòñÿ
				// óâåëè÷åíèåì ðàäèóñà ÷àñòèö.

				t_global = get_maximum_particle_time();

				for (long int i = 0; i < NP; ++i) {
					dt = t_global - particles[i].t;
					x = particles[i].x[0] + particles[i].vx * dt;

					if (DIFFUSE_RIGHT_WALL == true && x > L / 2 && particles[i].vx > 0)
						x = particles[i].x_max;

					if (x < min) min = x;
					if (x > max) max = x;
				}


				// Ðàññ÷èòûâàåì íà ñêîëüêî áëèçêî ÷àñòèöû íàõîäÿòñÿ ê ïåðâîé
				// è âòîðîé èäåàëüíûì ñòåíêàì

				min = min - particle_R;
				max = L - max - particle_R;

				// âûáèðàåì íàèìåíüøåå è ýòèõ ðàññòîÿíèé è ñîõðàíÿåì â deltaL
				deltaL = min;
				if (deltaL > max) deltaL = max;

				// ñæèìàåì íå âïðèòûê ê ÷àñòèöàì è íå ñëèøêîì áûñòðî
				deltaL = deltaL / 1.1;
				if (deltaL < 0.1e-12) deltaL = 0.01e-30; // íå äâèãàåì ñòåíêó åñëè ÷àñòèöû áëèçêî
				if (deltaL > delta_max) deltaL = delta_max;  // ñäâèãàåì ñòåíêó íå áîëüøå ÷åì íà delta_L

				dx = deltaL / 2.0;  // ðàññ÷èòûâàåì ñìåùåíèå äëÿ âñåõ ÷àñòèö


									// Åñëè ñëåäóùåå ñìåùåíèå ñòåíêè ñîæì¸ò ñèñòåìó áîëüøå
									// ÷åì òðåáóåòñÿ (ò.å. îòíîñèòåëüíàÿ ïëîòíîñòü ÷àñòèö â ñèñòåìå etta
									// áóäåò áîëüøå / ìåíüøå òîé, êîòîðàÿ áûëà óêàçàíà à àðãóìåíòàõ),
									// òî íåîáõîäèìî çàäàòü òî÷íîå çíà÷åíèå L, ÷òîáû ïëîòíîñòü ñîâïàëà
									// ñ òðåáóåìîé ïëîòíîñòüþ.

				if (L - deltaL < L_ideal) {
					deltaL = 0.0;
					L = L_ideal;
					dx = 0.0;
				}
			}
			else {

				//Åñëè ñèñòåìó íåîáõîäèìî ðàñøèðèòü, òî ïðîñòî áåð¸ì
				//íàèáîëüøèé ðàçðåø¸ííûé øàã äëÿ èçìåíåíèÿ êîîðäíàòû ñòåíêè
				//è ñäâèãàåì èäåàëüíóþ ñòåíêó èëè äâå ñòåíêè.

				deltaL = L - L_ideal;
				if (deltaL < -delta_max) deltaL = -delta_max; // øàã ðàñøèðåíèÿ ñèñòåìû
				dx = deltaL / 2.0;
			}

			if (type != 1) {

				// Äâèãàåì âñå ÷àñòèöû âëåâî, òàêèì îáðàçîì ïîäîäâèãàÿ òîëüêî ëåâóþ ñòåíêó.
				// ðàññòîÿíèå ÷àñòèö äî ïðàâîé ñòåíêè íå èçìåíèòñÿ, ò.ê. ïîñëå ïåðåìåùåíèÿ ÷àñòèö
				// âëåâî ìû ñäâèãàåì îáå ñòåíêè íà òî æå ðàññòîÿíèå - â èòîãå ìû ñäâèíåì âñå
				// ÷àñòèöû âëåâî íà deltaL/2.0.

				for (long int w = 0; w < NP; ++w) {
					for (int j = 0; j < 4; j++) {
						particles[w].x[j] -= dx;
					}
				}

				if (type == 0) {
					deltaL = dx;
				}
			}

			// èçìåíÿåì ñèñòåìó - ïðîèñõîäèò ìãíîâåííîå èçìåíåíèå êîîðäèíàò äâóõ ñòåíîê
			L -= deltaL;
		}
		else {
			double dy, dz, dt2, r, minR;
			max = -100.0;
			min = 1.1e+10;

			if (etta < compress_to_etta) {
				t_global = get_maximum_particle_time();
				minR = L;

				for (long int i = 0; i < NP; i++) {
					dt = t_global - particles[i].t;
					x = particles[i].x[0] + particles[i].vx * dt;

					if (x < min) min = x;
					if (x > max) max = x;

					for (long int j = 0; j < i; j++) {
						dt2 = t_global - particles[j].t;
						dx = particles[i].x[0] + particles[i].vx * dt - particles[j].x[0] - particles[j].vx * dt2;
						dy = particles[i].y[0] + particles[i].vy * dt - particles[j].y[0] - particles[j].vy * dt2;
						dz = particles[i].z[0] + particles[i].vz * dt - particles[j].z[0] - particles[j].vz * dt2;

						r = sqrt(dx*dx + dy*dy + dz*dz);

						if (r < minR) minR = r;
					}
				}

				minR = (minR - particle_Rx2) / 4.0;
				if (minR > delta_max) {
					minR = delta_max;
				}

				// Ñëåäèì ÷òîáû íå ïåðåñå÷ü ñòåíêè
				min = min - particle_R;
				max = L - particle_R - max;
				if (min > max) {
					min = max;
				}
				if (minR > min) {
					minR = min / 3.0;
				}

				if (minR < 1.0e-10) minR = 0.0;

				if (R_ideal > particle_R + minR) {
					change_R(particle_R + minR);
				}
				else {
					change_R(R_ideal);
				}
			}
			else {
				if (R_ideal < particle_R - delta_max / 2.0) {
					change_R(particle_R - delta_max / 2.0);
				}
				else {
					change_R(R_ideal);
				}
			}
		}

		if (temp_save == true) {
			save(temp_save_file);
		}

		//Ñîõðàíÿåì ñîñòîÿíèå ñèñòåìû è ñíîâà çàãðóæàåì åãî ïåðåñ÷èòàâ íîâûå
		//ïàðàìåòðû è ïðîâåäÿ íåîáõîäèìóþ èíèöèàëèçàöèþ
		save("tmp");
		load_seed("tmp");

		// Ïîñëå êàæäîãî ñìåùåíèÿ ñòåíêè äåëàåì óêàçàííîå êîëè÷åñòâî ñîóäàðåíèé
		// íà ÷àñòèöó â ñèñòåìå
		for (long int i = 0; i < KK1; i++) {
			step();
		}

		compression_steps_done += 1;
		if (compression_steps_done % KK2 == 0) {

			for (long int i = 0; i < KK3; i++) {
				step();

				// î÷èñòèòü ñòðîêó
				printf("\r%76c\r", ' ');

				printf("Relaxation: %ld / %ld, current etta=%.15le ", i, KK3, etta);
			}
		}

		// ðàññ÷èòûâàåì ïëîòíîñòü ïîñëå ñæàòèÿ
		etta = (4.0 * PI * particle_R3 * NP) / (3.0 * A * A * (L - particle_Rx2));

		// î÷èñòèòü ñòðîêó
		printf("\r%76c\r", ' ');

		printf("etta = %.15le, should be equal to %.15le ", etta, compress_to_etta);
	}

	printf("\n INFO: System density was sucessfully changed to %.15le \n", etta);
}

void measure_pressure(long int steps) {
	dV_left_wall = 0.0;
	dV_right_wall = 0.0;
	N_collissions_left_wall = 0.0;
	N_collissions_right_wall = 0.0;
	dP_time = 0.0;

	printf("\n Measure pressure on two ideal walls...\n");

	for (long int w = 1; w <= steps; w++) {
		// î÷èñòèòü ñòðîêó
		printf("\r%76c\r", ' ');
		printf("%ld / %ld", w, steps);
		step();
	}

	double m = 1.0;
	double P_left_wall = (2.0 * m * dV_left_wall) / (dP_time * A * A);
	double P_right_wall = (2.0 * m * dV_right_wall) / (dP_time * A * A);
	double N_collissions = (N_collissions_left_wall + N_collissions_right_wall) / 2.0;
	double accuracy = 1.0 / sqrtf(N_collissions);

	printf(" P_left_wall = %.15le \n P_right_wall = %.15le\n", P_left_wall, P_right_wall);
	printf("%ld collissions, accuracy = %f %% \n\n", long int(N_collissions), accuracy);

}

/*
Ôóíêöèÿ ïðîâåäåíèÿ ýêñïåðèìåíòà ïî óêàçàííîìó â
òåêñòîâîì ôàéëå îïèñàíèþ.

Àðãóìåíòû:
file_name - èìÿ ôàéëà ñ îïèñàíèåì øàãîâ ýêñïåðèìåíòà äëÿ ïðîãðàììû.
*/
void init(std::string file_name) {
	using namespace std;
	clock_t start, end, result;
	char command[255], parameter[255], tmp_str[255], parameter2[255];
	long i, steps;
	ifstream command_file(file_name.c_str());

	while (!command_file.eof()) {
		command_file >> command;
		string str_command = command;

		printf("\n\n<==========================>\n");
		srand(time(NULL));  // ïåðåìåøèâàåì ãåíåðàòîð ñëó÷àéíûõ ÷èñåë

							// Åñëè íåîáõîäèìî ñîçäàòü íîâûé ïîñåâ
		if (str_command.compare("new") == 0) {
			double etta;   // ñðåäíÿÿ ïëîòíîñòü ñèñòåìû
			long int particles_count;

			command_file >> particles_count;
			command_file >> etta;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			new_seed(particles_count, etta);

			print_system_parameters();
		}
		// Åñëè íåîáõîäèìî çàãðóçèòü ñîñòîÿíèå ñèñòåìû èç ôàéëà
		if (str_command.compare("load") == 0) {
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			load_seed(parameter);

			printf("\n System was successfully loaded from file '%s' \n", parameter);

			print_system_parameters();
		}
		// Åñëè íåîáõîäèìî ðàññ÷èòàòü äèíàìèêó ñèñòåìû â òå÷åíèè
		// óêàçàííîãî êîëè÷åñòâà ñîóäàðåíèé
		if (str_command.compare("step") == 0) {
			command_file >> steps;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			printf(" INFO: Step start. \n");

			start = clock();

			long int f[100];
			for (i = 0; i < 100; i++) f[i] = 0;
			double t_max;
			double dt, x;

			for (i = 0; i < steps; ++i) {
				step();

				// î÷èñòèòü ñòðîêó
				printf("\r%76c\r", ' ');

				printf("%ld / %ld", long(i + 1), steps);

				// åñëè ìû ïèøåì ôàéë ñ èñòîðèåé ñîáûòèé,
				// òî î÷èùàòü åãî êàæäûå 1000 ñîóäàðåíèé íà ÷àñòèöó,
				// ÷òîáû îí íå ïðåâûøàë äîïóñòèìûõ ðàçìåðîâ.
				/*
				if (i % 100000 == 0) {
				FILE *history_file = fopen("history.txt", "w+");
				fclose(history_file);

				check_particles();
				}
				*/

				if (DIFFUSE_RIGHT_WALL == true) {
					t_max = get_maximum_particle_time();

					for (long int j = 0; j < NP; j++) {
						dt = t_max - particles[j].t;
						x = particles[j].x[0] + particles[j].vx*dt;
						if (x > L - particle_R - x_diffuse) {
							f[long int(x * 100) / 100] += 1;
						}
					}

					if (i % 100 == 0) {
						long int mean = 0;
						for (int j = 0; j < 100; j++) {
							mean += f[j];
						}
						mean /= 100;

						for (long int j = 0; j < NP; j++) {
							dt = t_max - particles[j].t;
							x = particles[j].x[0] + particles[j].vx*dt;
							if (x > L - particle_R - x_diffuse) {
								if (f[long int(x * 100) / 100] > mean && particles[j].vx < 0.0) {
									particles[j].x_max = x;
								}
								else if (f[long int(x * 100) / 100] < mean && particles[j].vx < 0.0) {
									particles[j].x_max = x + (L - particle_R - x) / 2.0;
								}
							}
						}

						for (int j = 0; j < 100; j++) f[j] = 0;
					}
				}

			}

			end = clock();
			result = end - start;

			printf("\n INFO: finished %ld collisions per particle \n", steps);
			printf("Total Time = %f seconds. \n", double(result / CLOCKS_PER_SEC));

			/*
			Ðàññ÷èòûâàåì ìîìåíò èìïóëüñà ñèñòåìû Lx è âûâîäèì åãî íà ýêðàí.
			*/
			double Lx = 0.0;
			for (i = 0; i < NP; i++) {
				Lx += particles[i].y[0] * particles[i].vz - particles[i].z[0] * particles[i].vy;
			}
			printf("\n Lx = %.15le \n", Lx);
			check_particles();
		}
		// åñëè íåîáõîäèìî ñîáðàòü äàííûå ïî ïðîôèëþ ïëîòíîñòè â ñèñòåìå
		if (str_command.compare("image") == 0) {
			command_file >> steps; // ÷èñëî ñîóäàðåíèé íà îäíó ÷àñòèöó çà âðåìÿ èçìåðåíèé
			command_file >> i; // òî÷íîñòü. ÷èñëî òî÷åê ãðàôèêà íà îäèí ðàäèóñ ÷àñòèöû
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			image(steps, i, parameter);
		}
		// åñëè íåîáõîäèìî ñîáðàòü äàííûå ïî "ðàçðåçó" â ñèñòåìå
		if (str_command.compare("profile") == 0) {
			int dots_for_each_particle;
			double x1, x2;
			command_file >> x1;
			command_file >> x2;
			command_file >> dots_for_each_particle;
			command_file >> steps;
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			profile(x1, x2, dots_for_each_particle, steps, parameter);
		}
		if (str_command.compare("profile_medium") == 0) {
			int dots_for_each_particle;
			double x1, x2;
			command_file >> x1;
			command_file >> x2;
			command_file >> dots_for_each_particle;
			command_file >> steps;
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			profile_medium(x1, x2, dots_for_each_particle, steps, parameter);
		}
		if (str_command.compare("f1") == 0) {
			double x1, x2;
			int n1, n2, n3, n4;
			command_file >> x1;
			command_file >> x2;
			command_file >> n1;
			command_file >> n2;
			command_file >> n3;
			command_file >> n4;
			command_file >> parameter;

			printf("\n F1 started \n");

			F1(x1, x2, n1, n2, n3, n4, parameter);
		}
		// åñëè íåîáõîäèìî ñîõðàíèòü ñîñòîÿíèå ñèñòåìû
		if (str_command.compare("save") == 0) {
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			save(parameter);
			printf("\n INFO: coordinates of particles saved to '%s' \n", parameter);
		}
		if (str_command.compare("function_g") == 0) {
			double x1, x2;
			command_file >> steps;
			command_file >> x1;
			command_file >> x2;
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			printf("\n Function g started... %ld collissions \n", steps);
			function_g(steps, x1, x2, parameter);
			printf("\n Function g finished \n");
		}
		if (str_command.compare("count_of_nearest_particles") == 0) {
			double x1, x2;
			command_file >> steps;
			command_file >> x1;
			command_file >> x2;
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

													   //count_of_nearest_particles(steps, x1, x2, parameter);
			printf("\n Function count_of_nearest_particles finished \n");
		}
		if (str_command.compare("measure_pressure") == 0) {
			command_file >> steps;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			measure_pressure(steps);
		}
		if (str_command.compare("enable_autosave_for_compress") == 0) {
			temp_save = true;
			command_file >> temp_save_file;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè
		}
		if (str_command.compare("disable_autosave_for_compress") == 0) {
			temp_save = false;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè
		}
		// åñëè íåîáõîäèìî èçìåíèòü ïëîòíîñòü ñèñòåìû
		if (str_command.compare("compress") == 0) {
			command_file.getline(parameter, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè
			printf("\n Sorry, 'compress' command was deprecated, we need to use \n");
			printf("'compress_two_walls'(short form: 'compresst') instead. You can also use ");
			printf("'compress_left_wall'(short form: 'compressl') \n or 'compress_right_wall'");
			printf("(short form: 'compressr') commands to change system density.\n");
			exit(1);
		}
		// åñëè íåîáõîäèìî èçìåíèòü ïëîòíîñòü ñèñòåìû
		// ñæèìàòü, ñäâèãàÿ äâå èäåàëüíûå ñòåíêè
		if ((str_command.compare("compress_two_walls") == 0) ||
			(str_command.compare("compresst") == 0)) {
			long int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // òðåáóåìàÿ ïëîòíîñòü
			command_file >> delta_L;  // ìèíèìàëüíîå äîïóñòèìîå çíà÷åíèå ïåðåìåùåíèÿ ñòåíêè
			command_file >> KK1;  // êîëè÷åñòâî ñîóäàðåíèé ïîñëå êàæäîãî øàãà ñæàòèÿ
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè
			compress(etta, delta_L, KK1, KK2, KK3, 2);

			print_system_parameters();
		}
		// åñëè íåîáõîäèìî èçìåíèòü ïëîòíîñòü ñèñòåìû
		// äâèãàòü òîëüêî "ëåâóþ" èäåàëüíóþ ñòåêó, x = 0
		if ((str_command.compare("compress_left_wall") == 0) ||
			(str_command.compare("compressl") == 0)) {
			long int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // òðåáóåìàÿ ïëîòíîñòü
			command_file >> delta_L;  // ìèíèìàëüíîå äîïóñòèìîå çíà÷åíèå èçìåíåíèÿ ïëîòíîñòè
			command_file >> KK1;  // êîëè÷åñòâî ñîóäàðåíèé ïîñëå êàæäîãî øàãà ñæàòèÿ
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè
			compress(etta, delta_L, KK1, KK2, KK3, 0);

			print_system_parameters();
		}
		// åñëè íåîáõîäèìî èçìåíèòü ïëîòíîñòü ñèñòåìû
		// äâèãàòü òîëüêî "ïðàâóþ" èäåàëüíóþ ñòåêó, x = L
		if ((str_command.compare("compress_right_wall") == 0) ||
			(str_command.compare("compressr") == 0)) {
			long int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // òðåáóåìàÿ ïëîòíîñòü
			command_file >> delta_L;  // ìèíèìàëüíîå äîïóñòèìîå çíà÷åíèå èçìåíåíèÿ ïëîòíîñòè
			command_file >> KK1;  // êîëè÷åñòâî ñîóäàðåíèé ïîñëå êàæäîãî øàãà ñæàòèÿ
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè
			compress(etta, delta_L, KK1, KK2, KK3, 1);

			print_system_parameters();
		}
		// åñëè íåîáõîäèìî èçìåíèòü ïëîòíîñòü ñèñòåìû
		// ñæèìàòü, èçìåíÿÿ ðàäèóñ ÷àñòèö
		if (str_command.compare("compress_R") == 0) {
			long int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // òðåáóåìàÿ ïëîòíîñòü
			command_file >> delta_L;  // ìèíèìàëüíîå äîïóñòèìîå çíà÷åíèå èçìåíåíèÿ ðàäèóñà
			command_file >> KK1;  // êîëè÷åñòâî ñîóäàðåíèé ïîñëå êàæäîãî øàãà ñæàòèÿ
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè
			compress(etta, delta_L, KK1, KK2, KK3, 3);

			print_system_parameters();
		}
		if (str_command.compare("diffuse_right_wall_on") == 0 ||
			str_command.compare("diff_on") == 0) {
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			DIFFUSE_RIGHT_WALL = true;

			save("tmp.txt");
			load_seed("tmp.txt");
		}
		if (str_command.compare("diffuse_right_wall_off") == 0 ||
			str_command.compare("diff_off") == 0) {
			command_file.getline(tmp_str, 255, '\n');  // çàâåðøèòü ñ÷èòûâàíèå ñòðîêè

			DIFFUSE_RIGHT_WALL = false;

			save("tmp.txt");
			load_seed("tmp.txt");
		}
		if (str_command.empty()) break;
	}
	FILE *result_flag = fopen("result", "w+");
	fclose(result_flag);
}

/*
Ñòàðòîâàÿ ôóíêöèÿ äëÿ ïðîãðàììû, ñ÷èòàòü ôàéë "program.txt"
è íà÷àòü âûïîëíÿòü ýêñïåðèìåíò, îïèñàííûé â ýòîì ôàéëå.
*/
int main()
{
	FILE *history_file = fopen("history.txt", "w+");
	fclose(history_file);

	//particles_for_check[0] = 3026;
	particles_for_check_count = 0;

	init("program.txt");
	return 0;
}