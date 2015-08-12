/*
 1. не трогать функции работы с очередью событий, не оптимизировать.
    они отлажены и должны быть именно такими.
*/

/*
    ИМПОРТИРУЕМЫЕ МОДУЛИ
 */
#include "stdafx.h"
#include <time.h>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <exception>

using namespace System;

/*
    ОБЪЯВЛЕНИЕ ПЕРЕМЕНЫХ
 */

// длина стороны кубической ячейки, на которые разбивается объём 
short K, K2;

// число частиц во всём объёме
int NP = 6976;

int COLL_COUNT = 0;

// delete this dirty string
#define PI 3.141592653589793238462

// параметры объёма, задаются в load()
double A, A2, dA, L, dL, global_E = 0.0;

// индекс последнего элемента в очереди событий.
int last;

// объект "событие"
typedef struct Event_ {
    double t;
    int im, jm;
} Event;

// очередь событий - оптимально 8192 элемента
Event time_queue[30000];

// объект "частица"
// x,y,z,vx,vy,vz,t - particle coordinats.
// dt - delta time for the next event of this particle.
// kv - velocity koeficient.
typedef struct particle_ {
    double x, y, z, vx, vy, vz, t, dt, kv;
    int x_box, y_box, z_box, ti, box_i, i_copy;
} particle;

// массив частиц
particle particles[30000];

// клетка. Объём системы разделён на множество клеток,
// каждая клетка содержит в себе несколько виртуальных частиц
typedef struct Box_ {
    double x1, y1, z1, x2, y2, z2;
    int particles[11];
    short end;
} Box;

// массив клеток для всего объёма
Box boxes_yz[16][16][47];


// массив с номерами частиц, для которых надо сохранять историю событий
int particles_for_check[100];
int particles_for_check_count = 0;

int iim, ijm;


void print_time_line() {
    for (int i = 0; i <= last; i++) {
        fprintf(stderr, "EVENT:");
        fprintf(stderr, "%d %d %d %.5le", i, time_queue[i].im, time_queue[i].jm, time_queue[i].t);
    }
}


int check_particles() {
    double E = 0.0;

    //printf(" \n CHECK SYSTEM \n");

    for (int i = 0; i < NP; i++) {
        particle p1 = particles[i];
        Box p1_box = boxes_yz[p1.y_box][p1.z_box][p1.x_box];

        E += p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;

        if ((p1.x_box > K2) || (p1.y_box > K) || (p1.z_box > K) ||
            (p1.x_box < 0) || (p1.y_box < 0) || (p1.z_box < 0)) {
            return 400;
        }

        if (((abs(p1.x) - L) > 1.0e-14) ||
			((abs(p1.y) - A) > 1.0e-14) ||
			((abs(p1.z) - A) > 1.0e-14)) {
            printf(" \n Particle %d, %.15le, %.15le, %.15le, A = %.15le, %.15le \n ", i, p1.x, p1.y, p1.z, A, p1.z - A);
            return 500;
        }

        if ((p1.x < p1_box.x1) || (p1.x > p1_box.x2) ||
            (p1.y < p1_box.y1) || (p1.y > p1_box.y2) ||
            (p1.z < p1_box.z1) || (p1.z > p1_box.z2)) {

            printf("Vilet za granicy %d \n", i);

            printf("Granizy:\n");
            printf("X : [%.15le ; %.15le]\n", p1_box.x1, p1_box.x2);
            printf("Y : [%.15le ; %.15le]\n", p1_box.y1, p1_box.y2);
            printf("Z : [%.15le ; %.15le]\n", p1_box.z1, p1_box.z2);
            printf("x, y, z: %.15le %.15le %.15le\n", p1.x, p1.y, p1.z);

            printf("p1.dt = %.15le, p.im = %d, p1.jm = %d \n", p1.t, time_queue[p1.ti].im, time_queue[p1.ti].jm);

            return 1;
        }

        if ((p1.i_copy > -1) && (particles[i+NP].i_copy == -1)) return 2;

        bool w = false;
        for (int ty = 0; ty <= boxes_yz[p1.y_box][p1.z_box][p1.x_box].end; ty++) {
            //printf(" %d ", boxes_yz[p1.y_box][p1.z_box][p1.x_box].particles[ty]);
            if (boxes_yz[p1.y_box][p1.z_box][p1.x_box].particles[ty] == i) w = true;
        }
        if (w == false) {
            for (int t = 0; t <= boxes_yz[p1.y_box][p1.z_box][p1.x_box].end; t++) {
                printf(" %d ", boxes_yz[p1.y_box][p1.z_box][p1.x_box].particles[t]);
            }

            return 1010;
        }
        if (time_queue[p1.ti].im != i && time_queue[p1.ti].jm != i) {
            printf("\n >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ");
            printf("\n i = %d ; im = %d ; jm = %d ; ti = %d ", i, time_queue[p1.ti].im, time_queue[p1.ti].jm, p1.ti);
            return 3;
        }
    }

    if (abs(E - global_E) > 0.1e-8) {
        printf("\nENERGY was changed: \n E_seed = %.15le \n E_now= %.15le \n", global_E, E);
        return 1000;
    }

	double Mx = 0.0;
	for (int i = 0; i < NP; i++) {
		Mx += (particles[i].y)*particles[i].vz - (particles[i].z)*particles[i].vy;
	}
	if (abs(Mx) > 1.0e-14) {
		printf(" \n Mx = %.15le \n", Mx);
		return 2501;
	}
///////////////////// ?????????????????????????????????????????????????????????????????????????????????????????????
    for (int i = 0; i < last; i++) {
        if ((time_queue[i].im >= NP) && (particles[time_queue[i].im-NP].i_copy == -1)) {
            printf("\n Time tree event # %d ", i);
            printf("\n im, jm = %d %d \n", time_queue[i - 2].im, time_queue[i - 2].jm);
            printf("\n im, jm = %d %d \n", time_queue[i-1].im, time_queue[i-1].jm);
            printf("\n im, jm = %d %d \n", time_queue[i].im, time_queue[i].jm);
            printf("\n im, jm = %d %d \n", time_queue[i+1].im, time_queue[i+1].jm);
            printf("\n particle 34 event # %d \n", particles[34].ti);
            printf("\n last = %d \n", last);
            return time_queue[i].im;
        }
        if ((time_queue[i].jm >= NP) && (particles[time_queue[i].jm-NP].i_copy == -1)) {
            return time_queue[i].im;
        }
    }

    return 0;
}


/*
 функция подъема элемента по очереди событий к началу очереди
 параметры:
 i - позиция, на которой находится элемент в данный момент
 t - время до наступления данного события
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
        } else return i;
    }
    return 1;
}

/*
 функция добавления события в очередь событий
 параметры:
 e - новое событие (см. подробнее структуру типа Event)
 */
void add_event(int &i, int &j) {
    double t = particles[i].t + particles[i].dt;

    particles[i].ti = get_up(last, t);

    if (j >= 0) {
        particles[j].dt = particles[i].dt;
        particles[j].ti = particles[i].ti;
    }

    time_queue[ particles[i].ti ].im = i;
    time_queue[ particles[i].ti ].jm = j;
    time_queue[ particles[i].ti ].t = t;

    last++;
}

/*
 функция удаления события из очереди событий
 параметры:
 i - позиция удаляемого элемента в очереди событий
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
        } else {
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
    int k = particles[i].ti;

	///////////////
	/*
	for (int t = 0; t < particles_for_check_count; t++) {
		if (particles_for_check[t] == time_queue[k].im || particles_for_check[t] == time_queue[k].jm) {
			FILE *save_file = fopen("history.txt", "a");
			fprintf(save_file, "\n  CLEAR EVENT: %d %d \n", time_queue[k].im, time_queue[k].jm);
			fclose(save_file);
		}
	}
	*/
	///////////////

    if (k > 0) {
		int e = -100;
		int kim = time_queue[k].im;
		int kjm = time_queue[k].jm;

        if (time_queue[k].im == i) {
            if (time_queue[k].jm >= 0) {
                //time_queue[k].im = time_queue[k].jm;
                //time_queue[k].jm = -100;

				double dt = particles[kim].t - particles[kjm].t;  // разница в текущем времени

				particles[kjm].x += particles[kjm].vx*dt;
				particles[kjm].y += particles[kjm].vy*dt;
				particles[kjm].z += particles[kjm].vz*dt;
				particles[kjm].t = particles[kim].t;
				particles[kjm].dt = (particles[kjm].dt - dt) / 1.1;

				delete_event(k);
				add_event(kjm, e);
            }
            else delete_event(k);
        }
        else if (time_queue[k].jm == i) {
			//time_queue[k].jm = -100;

			double dt = particles[kjm].t - particles[kim].t;  // разница в текущем времени

			particles[kim].x += particles[kim].vx*dt;
			particles[kim].y += particles[kim].vy*dt;
			particles[kim].z += particles[kim].vz*dt;
			particles[kim].t = particles[kjm].t;
			particles[kim].dt = (particles[kim].dt - dt) / 1.1;

			delete_event(k);
			add_event(kim, e);
        }
        particles[i].ti = 0;
    }
}


/*
 функция расчёта ближайшего события для частицы
 параметры:
 i - номер частицы
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
    } else {
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
    } else {
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
    } else {
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
                if (q == K+1) {
                    continue;
                }
                if (w == -1) {
                    continue;
                }
                if (w == K+1) {
                    continue;
                }

                for (s = 0; s <= boxes_yz[q][w][r].end; ++s) {
                    n = boxes_yz[q][w][r].particles[s];

                    if (n == p1.i_copy) continue;

                    temp = p1.t - particles[n].t;

                    dvx = particles[n].vx - p1.vx;
                    dvy = particles[n].vy - p1.vy;
                    dvz = particles[n].vz - p1.vz;
                    dx = particles[n].x + particles[n].vx * temp - p1.x;
                    dy = particles[n].y + particles[n].vy * temp - p1.y;
                    dz = particles[n].z + particles[n].vz * temp - p1.z;

                    bij = dx * dvx + dy * dvy + dz*dvz;
                    if (bij < 0.0) {
                        dv = dvx * dvx + dvy * dvy + dvz*dvz;
                        d = bij * bij + dv * (4.0 - dx * dx - dy * dy - dz * dz);

                        if (d > 0.0) {
                            dt = -(sqrt(d) + bij) / dv;

							/*
							if (dt < -1.0e-12) {
								double d = dx*dx + dy*dy + dz*dz;
								if (i < NP && n < NP && d < 4.0 - 1.0e-14) {
									printf("\n DISTANCE = %.15le", d);
									printf("\n particle %d: vx = %.15le, vy = %.15le, vz = %.15le \n", i, particles[i].vx, particles[i].vy, particles[i].vz);
									printf("particle %d: vx = %.15le, vy = %.15le, vz = %.15le \n", n, particles[n].vx, particles[n].vy, particles[n].vz);
									printf("\n i = %d, n = %d, i.i_copy = %d, n.i_copy = %d \n", i, n, particles[i].i_copy, particles[n].i_copy);
									printf("11");
								}
							}
							if (dt < -1.0e-15 && dt > -1.0e-12) {
								printf("\n previous event: iim = %d, ijm = %d \n", iim, ijm);
								printf("\n i = %d, n = %d, dt = %.15le, distance = %.15le \n", i, n, dt, dx*dx + dy*dy + dz*dz);
							}
							*/
							/*
							 Сценарии, при котором возможны отрицательные времена:
							  1. Образ, вставленный в систему для мгновенного соударения может проникать
							     в другие частицы, после мгновенного соударения с частицей из объема этот
								 образ будет уничтожен.
							  2. Частица n1 сталкивается с частицей n2, мы создаем образ для частицы n1,
							     для которого не находится места и мы сталкиваем его с другой частицей n3,
							     в результате чего скорость частицы n1 снова меняется и время соударения
							     частиц n1 и n2 может быть отрицательным (-1*10-14) из за погрешности
							     в расчете координат в 15ом знаке.
							*/
							if (dt > -1.0e-12 && dt < 1.0e-15) dt = 0.0;

                            temp += dt;

							/*
                            for (int t = 0; t < particles_for_check_count; t++) {
                                if (i == particles_for_check[t] || n == particles_for_check[t]) {
                                    FILE *save_file = fopen("history.txt", "a");
                                    fprintf(save_file, "\n retime tempory:: i = %d, n = %d, dt = %.15le \n", i, n, dt);
                                    fprintf(save_file, "particle %d x = %.15le, p1.x = %.15le", i, particles[i].x, p1.x);
                                    fclose(save_file);
                                }
                            }
							*/

                            if ((dt < dt_min) && (dt >= 0.0) &&
                            ((temp < particles[n].dt) || ((abs(temp - particles[n].dt) < 0.1E-15) &&
                            (time_queue[particles[n].ti].im == n) &&
                            (time_queue[particles[n].ti].jm == -100)))) {

                                ////////////////////
								/*
                                for (int t = 0; t < particles_for_check_count; t++) {
                                    if (i == particles_for_check[t] || n == particles_for_check[t]) {
                                        FILE *save_file = fopen("history.txt", "a");
                                        fprintf(save_file, "\n ****** dt_min = %.15le \n ", dt);
                                        fclose(save_file);
                                    }
                                }
								*/
                                ////////////////////

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

		//////////////////
		/*
		for (int t = 0; t < particles_for_check_count; t++) {
			if (particles_for_check[t] == i || particles_for_check[t] == jm || particles_for_check[t] == particles[i].i_copy) {
				FILE *save_file = fopen("history.txt", "a");
				fprintf(save_file, "\n  New event: im = %d, jm = %d, dt = %.15le, delta T between particles = %.15le \n", i, jm, dt_min, dt);
				fclose(save_file);
			}
		}
		*/
		//////////////////

        particles[jm].t = p1.t;
        particles[jm].dt = dt_min;
        particles[jm].x += particles[jm].vx * dt;
        particles[jm].y += particles[jm].vy * dt;
        particles[jm].z += particles[jm].vz * dt;

        clear_particle_events(jm);
    }

    particles[i].dt = dt_min;
    add_event(i, jm);

    //////////////////
	/*
    for (int t = 0; t < particles_for_check_count; t++) {
        if (particles_for_check[t] == i || particles_for_check[t] == jm || particles_for_check[t] == particles[i].i_copy) {
			FILE *save_file = fopen("history.txt", "a");
            fprintf(save_file, "\n    retime result: %d %d %.15le \n", i, jm, dt_min);
			fclose(save_file);
        }
    }
	*/
    //////////////////

    if (dt_min < -1.0e-11) {
		printf("\n p1.x = %.15le, p1.y = %.15le, p1.z = %.15le \n", p1.x, p1.y, p1.z);
		printf("\n p1.x_box = %d, p1.y_box = %d, p1.z_box = %d", p1.x_box, p1.z_box, p1.y_box);
		printf("\n im = %d, jm = %d, dt = %.15le, A = %.15le", i, jm, dt_min, A);
        printf("QQQQQQQQQQQQQQ");
    }
}


// Данная функция позволяет найти свободное место для нового образа
// по оси Х. Если место не найдено, функция сталкивает новый образ
// с одной из частиц, мешающих его поставить
void find_place_for_particle(int &i) {

    double x, y, z, dy, dz, dx, d, dt, dvx, dvy, dvz, bij, dv, x_min, x_max = 0.0;
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

    for (int x_box = 0; x_box <= K2; ++x_box)
        for (int y_box = particles[i].y_box - 1; y_box < particles[i].y_box + 2; ++y_box) {
            for (int z_box = particles[i].z_box - 1; z_box < particles[i].z_box + 2; ++z_box) {

                if (y_box == -1) {
                    continue;
                }
                if (y_box == K+1) {
                    continue;
                }
                if (z_box == -1) {
                    continue;
                }
                if (z_box == K+1) {
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
        for (int r2 = r1+1; r2 <= spaces; r2++)
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
        // if two spaces overlaps with each other or if we have no enough free space 
        // between two areas, we need to combine them
        //abs(no_free_space_max[r1] - no_free_space_min[r1 + 1]) < 1.0e-12)
        if (no_free_space_max[r1] > no_free_space_min[r1+1])
        {
            if (no_free_space_max[r1] < no_free_space_max[r1+1])
                no_free_space_max[r1] = no_free_space_max[r1+1];
            for (int r2 = r1+1; r2 < spaces; r2++)
            {
                no_free_space_min[r2] = no_free_space_min[r2+1];
                no_free_space_max[r2] = no_free_space_max[r2+1];
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
            particle p1 = particles[i];
            particle p2;
            double fa, fb, fc, fD, sD, t1, t2;

            for (int j = 0; j < NP; j++) {

                if (j == i - NP) continue;

                p2 = particles[j];

                dt = p1.t - p2.t;

                // if we have small time differences for these two particles
                // we shouldn't use it for collission because other virtual particle
                // can use the same particle in the same time.
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

                            return;
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

                                return;
                            }
                        }
                    }
                }
            }

			printf("\n A = %.15le, L = %.15le\n ", A, L);
			printf("\n %d particle: x = %.15le, y = %.15le, z = %.15le \n", i, particles[i].x, particles[i].y, particles[i].z);
			printf("\n %d particle: vx = %.15le, vy = %.15le, vz = %.15le \n", i, particles[i].vx, particles[i].vy, particles[i].vz);
            printf("\n");
        }
    }

    x_max = -L;
    x_min = L;
    double length = 2*L, r;

    // select minimal free space and insert new particle in this space
    for (int m = 1; m < spaces; ++m) {
        // free space distance:
        r = abs(no_free_space_max[m] - no_free_space_min[m+1]);
        if (r < length) {
            length = r;
            x_min = no_free_space_max[m];
            x_max = no_free_space_min[m+1];
        }
    }

    //FILE *save_file2 = fopen("history.txt", "a");
    //fprintf(save_file2, "\n x_min: %.15le, x_max: %.15le \n", x_min, x_max);
    //fclose(save_file2);

    particles[i].x = x_min + (x_max - x_min)/2.0;
}


void destroy_virt_particle(int &i) {
    if (i >= NP) {
        printf("\n %d particle \n", i);
        printf(" i_copy = %d   particle icopy = %d \n", particles[i].i_copy, particles[i-NP].i_copy);
        printf(" x, y, z: %.5le %.5le %.5le", particles[i].x, particles[i].y, particles[i].z);

        throw "AAAAAAAAAAA";
    }
    if (particles[i].i_copy == -1) return;

    int new_i = i + NP;
    int x_box, y_box, z_box, box_i;

    x_box = particles[new_i].x_box;
    y_box = particles[new_i].y_box;
    z_box = particles[new_i].z_box;
    box_i = particles[new_i].box_i;
    Box p_box = boxes_yz[y_box][z_box][x_box];

    //////////
	/*
    for (int t = 0; t < particles_for_check_count; t++) {
        if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
            FILE *save_file = fopen("history.txt", "a");
            fprintf(save_file, "particle %d should be destroyed... \n", new_i);
            fclose(save_file);
        }
    }
	*/
    //////////

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

    //printf("particle %d was destroyed...", new_i);

    //////////
	/*
    for (int t = 0; t < particles_for_check_count; t++) {
        if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
            FILE *save_file = fopen("history.txt", "a");
            fprintf(save_file, "particle %d was destroyed... \n", new_i);
            fclose(save_file);
        }
    }
	*/
    //////////
}



void change_with_virt_particles(int &im, int &jm) {
    int y_box, z_box, box_i;
    int f = im + NP;

    //printf("\n im f = %d %d %d \n", im, f, particles[f].ti);
	//printf("\n p.i_copy = %d \n", particles[im].i_copy);

    double dt = particles[im].t - particles[f].t;
    particles[im].x = particles[f].x + dt*particles[f].vx;
    particles[im].y = particles[f].y + dt*particles[f].vy;
    particles[im].z = particles[f].z + dt*particles[f].vz;

    // записываем её в новую ячейку
    y_box = particles[f].y_box;
    z_box = particles[f].z_box;

    if (y_box == 0) y_box = 1;
    if (y_box == K) y_box = K - 1;
    if (z_box == 0) z_box = 1;
    if (z_box == K) z_box = K - 1;

    particles[im].x_box = particles[f].x_box;
    particles[im].y_box = y_box;
    particles[im].z_box = z_box;

    // при обмене виртуального образа на частицу ставим частицу точно на границу ячейки,
    // чтобы избегать накопления ошибок
    if (particles[im].y > A) particles[im].y = boxes_yz[y_box][z_box][particles[f].x_box].y2;
    else if (particles[im].y < -A) particles[im].y = boxes_yz[y_box][z_box][particles[f].x_box].y1;
    if (particles[im].z > A) particles[im].z = boxes_yz[y_box][z_box][particles[f].x_box].z2;
    else if (particles[im].z < -A) particles[im].z = boxes_yz[y_box][z_box][particles[f].x_box].z1;

    //printf("... Change with virt particle destroy ...");
	/*
    for (int t = 0; t < particles_for_check_count; t++) {
        if (im == particles_for_check[t]) {
            FILE *save_file = fopen("history.txt", "a");
            fprintf(save_file, "particle %d Change with virt particle destroy... \n", im);
            fclose(save_file);
        }
    }
	*/

    destroy_virt_particle(im);

    //clear_particle_events(im);
}


void create_virt_particle(int &i) {
    double dt, dt_min, t1, t2, t11, t22, t_min, y, z, dy, dz;
    double x2, y2, z2, kv = 1.0;
    double t01, t02;
    short x_box, y_box, z_box;
    int new_i = i + NP;

    ////////////////////////////
	/*
    printf("... Create virt particle destroy ...");
    for (int t = 0; t < particles_for_check_count; t++) {
        if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
            FILE *save_file = fopen("history.txt", "a");
            fprintf(save_file, "particle %d Create virt particle destroy... \n", i);
            fclose(save_file);
        }
    }
	*/
    ////////////////////////////

    destroy_virt_particle(i);

    ////////////////////////////
	/*
    printf("Trying to create virt particle %d...", new_i);
    for (int t = 0; t < particles_for_check_count; t++) {
        if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
            FILE *save_file = fopen("history.txt", "a");
            fprintf(save_file, "Trying to create virt particle %d... \n", new_i);
            fclose(save_file);
        }
    }
	*/
    ////////////////////////////

    if (((particles[i].y <= 1.0 - A) && (particles[i].vy < 0.0)) ||
        ((particles[i].y >= A - 1.0) && (particles[i].vy > 0.0)) ||
        ((particles[i].z <= 1.0 - A) && (particles[i].vz < 0.0)) ||
        ((particles[i].z >= A - 1.0) && (particles[i].vz > 0.0))) {

        ////////////////////////////
		/*
        printf("Starting to create virt particle %d...", new_i);
        for (int t = 0; t < particles_for_check_count; t++) {
            if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
                FILE *save_file = fopen("history.txt", "a");
                fprintf(save_file, "Starting to create virt particle %d... \n", new_i);
                fclose(save_file);
            }
        }
		*/
        ////////////////////////////

        dt_min = 1.0e+20;

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

        if ((particles[i].y*particles[i].vy > 0) && (dy <= 1.0)) {
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

                //printf("\n %d Rozdenie 1go roda \n", i);

                ////////////////////////////
				/*
                for (int t = 0; t < particles_for_check_count; t++) {
                    if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
                        FILE *save_file = fopen("history.txt", "a");
                        fprintf(save_file, "\n %d Rozdenie 1go roda \n", i);
                        fclose(save_file);
                    }
                }
				*/
                ////////////////////////////
            }
            else {
                t_min = t02;

                particles[new_i].x = particles[i].x - particles[i].vx*t02 - particles[new_i].vx*dt;
                particles[new_i].y = particles[i].y - particles[i].vy*t02 - particles[new_i].vy*dt;
                particles[new_i].z = particles[i].z - particles[i].vz*t02 - particles[new_i].vz*dt;

                ////////////////////////////
				/*
                for (int t = 0; t < particles_for_check_count; t++) {
                    if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
                        FILE *save_file = fopen("history.txt", "a");
                        fprintf(save_file, "\n %d Rozdenie 2go roda \n", i);
                        fclose(save_file);
                    }
                }
				*/
                ////////////////////////////
            }

            dt_min = dt;
        }
        if ((particles[i].z*particles[i].vz > 0) && (dz <= 1.0)) {
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

                    ////////////////////////////
					/*
                    for (int t = 0; t < particles_for_check_count; t++) {
                        if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
                            FILE *save_file = fopen("history.txt", "a");
                            fprintf(save_file, "\n %d Rozdenie 3go roda \n", i);
                            fclose(save_file);
                        }
                    }
					*/
                    ////////////////////////////
                }
                else {
                    t_min = t02;

                    particles[new_i].x = particles[i].x - particles[i].vx*t02 - particles[new_i].vx*dt;
                    particles[new_i].y = particles[i].y - particles[i].vy*t02 - particles[new_i].vy*dt;
                    particles[new_i].z = particles[i].z - particles[i].vz*t02 - particles[new_i].vz*dt;

                    ////////////////////////////
					/*
                    for (int t = 0; t < particles_for_check_count; t++) {
                        if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
                            FILE *save_file = fopen("history.txt", "a");
                            fprintf(save_file, "\n %d Rozdenie 4go roda \n", i);
                            fclose(save_file);
                        }
                    }
					*/
                    ////////////////////////////
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
            int k = x / (2.0 * LL);
            x = x - 2.0 * k * LL;
            int x_wall = 1;
            if (k % 2 != 0)
                x_wall = -1;
            else
                particles[new_i].vx = -particles[new_i].vx;
            particles[new_i].x = x_wall * (f * LL - x);
        }

        x_box = short((L + dL + particles[new_i].x) / dL);
        y_box = short((A + dA + particles[new_i].y) / dA);
        z_box = short((A + dA + particles[new_i].z) / dA);
        
        if (y_box < 0) y_box = 0;
        if (y_box > K) y_box = K;
        if (z_box < 0) z_box = 0;
        if (z_box > K) z_box = K;

        particles[new_i].x_box = x_box;
        particles[new_i].y_box = y_box;
        particles[new_i].z_box = z_box;

        particles[new_i].t = particles[i].t;

        // проверка того что новый образ не пересекается с уже существующими частицами
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
                            //printf("\n SEARCH!!! %d %d \n", new_i, n);

                            ////////////////////////////
							/*
                            for (int t = 0; t < particles_for_check_count; t++) {
                                if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
                                    FILE *save_file = fopen("history.txt", "a");
                                    fprintf(save_file, "\n SEARCH!!! %d %d \n", new_i, n);
                                    fclose(save_file);
                                }
                            }
							*/
                            ////////////////////////////
                        }
                    }
                }
            }
        }

        ////////////////////////////
		/*
        for (int t = 0; t < particles_for_check_count; t++) {
            if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
                FILE *save_file = fopen("history.txt", "a");
                fprintf(save_file, "\n Before find: x = %.15le, y = %.15le, z = %.15le \n", particles[new_i].x, particles[new_i].y, particles[new_i].z);
                fclose(save_file);
            }
        }
		*/
        ////////////////////////////

        if (search == true) {
            find_place_for_particle(new_i);

            /* update information about cell for new virtual particle */
            x_box = short((L + dL + particles[new_i].x) / dL);
            y_box = short((A + dA + particles[new_i].y) / dA);
            z_box = short((A + dA + particles[new_i].z) / dA);

            if (y_box < 0) y_box = 0;
            if (y_box > K) y_box = K;
            if (z_box < 0) z_box = 0;
            if (z_box > K) z_box = K;

            particles[new_i].x_box = x_box;
            particles[new_i].y_box = y_box;
            particles[new_i].z_box = z_box;
        }

        ////////////////////////////
		/*
        for (int t = 0; t < particles_for_check_count; t++) {
            if (i == particles_for_check[t] || new_i == particles_for_check[t]) {
                FILE *save_file = fopen("history.txt", "a");
                fprintf(save_file, "\n After find: x = %.15le, y = %.15le, z = %.15le \n", particles[new_i].x, particles[new_i].y, particles[new_i].z);
                fclose(save_file);
            }
        }
		*/
        ////////////////////////////

        
        short end = particles[new_i].box_i = ++boxes_yz[y_box][z_box][x_box].end;
        boxes_yz[y_box][z_box][x_box].particles[end] = new_i;

        particles[new_i].i_copy = i;
        particles[i].i_copy = new_i;
    }
}



/*
 функция загрузки данных о системе из файла.
 загружаются параметры объёма, координаты и скорости всех частиц
 */
void load_seed(std::string file_name) {
    double a1, a2, a3, x, y, z;
    int i;
    FILE *loading_file = fopen(file_name.c_str(), "r");
    fscanf(loading_file, "%i\n", &i);
    NP = i;
    fscanf(loading_file, "%le\n", &a1);
    A = a1 / 2.0;
    A2 = a1;
    fscanf(loading_file, "%le\n", &a1);
    L = a1 / 2.0;
    dA = A2 / (K - 1);
    dL = L * 2.0 / (K2 - 1);

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

    for (int i = 0; i < NP*2; i++) {
        particles[i].x_box = 0;
        particles[i].y_box = 0;
        particles[i].z_box = 0;
        particles[i].box_i = 0;
        particles[i].i_copy = -1;
        particles[i].ti = 0;
        particles[i].dt = 0.0;
    }

	global_E = 0.0;
    short end, x_box, y_box, z_box;
    for (int i = 0; i < NP; i++) {
        fscanf(loading_file, "%le %le %le\n", &a1, &a2, &a3);
        particles[i].x = a1 - L;
        particles[i].y = a2 - A;
        particles[i].z = a3 - A;
        x_box = short((L + dL + particles[i].x) / dL);
        y_box = short((A + dA + particles[i].y) / dA);
        z_box = short((A + dA + particles[i].z) / dA);

        if (boxes_yz[y_box][z_box][x_box].x1 <= particles[i].x &&
                boxes_yz[y_box][z_box][x_box].x2 >= particles[i].x)
            particles[i].x_box = x_box;
        else {
            printf("Particle locates in incorrect place %d", i);
            throw "Particle locates in incorrect place";
        }

        if (boxes_yz[y_box][z_box][x_box].y1 <= particles[i].y &&
                boxes_yz[y_box][z_box][x_box].y2 >= particles[i].y)
            particles[i].y_box = y_box;
        else {
            printf("Particle locates in incorrect place %d", i);
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
            throw "Error: Too many particles in one box";
        }
        for (short t = 0; t < particles[i].box_i; ++t)
            for (short d = t + 1; d <= particles[i].box_i; ++d)
                if (boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].particles[t] == boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].particles[d])
                {
                    printf("Duplicated particles in one box");
                    throw "Duplicated particles in one box";
                }
    }
    fclose(loading_file);

    last = 1;
    time_queue[0].t = 0.0;
    time_queue[0].im = -1;
    for (int i = 1; i < 30000; ++i) time_queue[i].t = 1.0E+20;

    for (int i = 0; i < NP; ++i) retime(i);
	for (int i = 0; i < NP; ++i) {
		create_virt_particle(i);
		if (particles[i].i_copy > 0) retime(particles[i].i_copy);
	}

	double Mx = 0.0;
	for (int i = 0; i < NP; i++) {
		Mx += (particles[i].y + A)*particles[i].vz - (particles[i].z + A)*particles[i].vy;
	}

	printf("\n Mx = %.15le\n", Mx);
}

/*
 функция сохранения состояния системы в текстовый файл.
 перед сохранением производится синхронизация частиц по времени.
 параметры:
 file_name - имя файла, в который будет записана информация
 */
void save(std::string file_name) {
    double x, y, z;

    FILE *save_file = fopen(file_name.c_str(), "w+");
    fprintf(save_file, "%d\n", NP);
    fprintf(save_file, "%.15le\n", A * 2.0);
    fprintf(save_file, "%.15le\n", L * 2.0);
    for (int i = 0; i < NP; ++i) {
        x = L + particles[i].x - particles[i].vx * particles[i].t;
        y = A + particles[i].y - particles[i].vy * particles[i].t;
        z = A + particles[i].z - particles[i].vz * particles[i].t;
        fprintf(save_file, "%.15le %.15le %.15le\n", x, y, z);
        fprintf(save_file, "%.15le %.15le %.15le\n", particles[i].vx, particles[i].vy, particles[i].vz);

        if (-particles[i].t > particles[i].dt)
            printf("Error: Sync of system events failed\n");
    }
    fclose(save_file);
}


// процедура нового посева системы.
// получает в качестве параметров
// NN - число частиц в ребре
// etta - начальная плотность
void new_seed(int NN, double etta) {
    int KZ = 2;

    double axy, axz, v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z;

    double rb = 2.0 * sqrt(2.0); // расстояние между двумя слоями
    double rbp = sqrt(2.0);
    double A0 = rb * (NN - 1.0) + 2.0; // рассчет параметров начального объема
    double L0 = rb * (2.0 * NN - 1.0) + 2.0; //

    double mL0 = L0;

    double Betta = 0.0;

    double XC, YC, ZC; // координаты центра объема
    int NA; // количество посеянных частиц
    double ax, ay, az, vx, vy, vz;

    NP = 2.0 * NN * (NN * NN + (NN - 1.0) * (NN - 1.0)); //  ожидаемое число частиц.
    NP = NP + (2.0 * NN - 1.0) * 2.0 * (NN - 1.0) * NN; //

    double etta0 = (4.0 * PI * NP) / (3.0 * A0 * A0 * (L0 - 2.0)); // рассчет первоначальной плотности.

    double Alpha = etta0 / etta;

    double sk = L0 / (2.0 * Alpha * (L0 - 2.0));

    Betta = exp((1.0 / 3.0) * log(L0 / (2.0 * Alpha * (L0 - 2.0)) + sk));
    //Betta = Betta + exp( (1.0 / 3.0) * log( -L0 / (2.0 * Alpha * (L0-2.0)) - sk) );

    Betta = 1.0 / Betta;

    XC = L0 / 2.0; //
    YC = A0 / 2.0; //  вычисление центра объема.
    ZC = A0 / 2.0; //

    L0 = L0*KZ;
    L = L0 * Betta;  // множитель Betta - это коэффициент расширения
    A = A0 * Betta;  // на него умножаются все координаты и параметры объема

    NA = 0;

	srand(501);

    for (int i = 0; i < 2 * NN; i++)
        for (int j = 0; j < NN; j++)
            for (int k = 0; k < NN / 2; k++) {
                vx = double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)) - double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)); //  задаем скорость.
                vy = double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)) - double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)); //  она будет одинакова для нескольких частиц сразу(с небольшим отклонением)
                vz = double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)) - double(rand() % 1000 + 1) / (1000.0 + double(rand() % 100)); //  чтобы был равен нулю импульс и момент импульса системы.

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

                            particles[NA].x = ax; // назначаем координаты для новой молекулы
                            particles[NA].y = ay;
                            particles[NA].z = az;

                            if ((j == 0) && (fabs(ZC - az) < 0.1E-15L) && ((i / 2)*2 != i))
                                if (((i + 1) / 4)*4 != (i + 1)) {
                                    axy = vy;
                                    axz = vz;
                                    if (ii > 0) {
                                        vy = -vy;
                                        vz = -vz;
                                    }
                                    vy = ii*vy;
                                    vz = ii*vz;
                                } else {
                                    vy = axy;
                                    vz = axz;
                                    vy = ii*vy;
                                    vz = ii*vz;
                                }
                            if (i == 0) {
                                if ((j == 0) && (kk < 0))
                                    switch (k) {
                                        case 0: v0x = vx;
                                            v0y = vy;
                                            v0z = vz;
                                            break;
                                        case 1: v1x = vx;
                                            v1y = vy;
                                            v1z = vz;
                                            break;
                                        case 2: v2x = vx;
                                            v2y = vy;
                                            v2z = vz;
                                            break;
                                        case 3: v3x = vx;
                                            v3y = vy;
                                            v3z = vz;
                                            break;
                                    } else
                                    switch (k) {
                                        case 0: vx = -v0x;
                                            vy = -v0y;
                                            vz = -v0z;
                                            break;
                                        case 1: vx = -v1x;
                                            vy = -v1y;
                                            vz = -v1z;
                                            break;
                                        case 2: vx = -v2x;
                                            vy = -v2y;
                                            vz = -v2z;
                                            break;
                                        case 3: vx = -v3x;
                                            vy = -v3y;
                                            vz = -v3z;
                                            break;
                                    }

                                if ((k == 0) && ((j / 2)*2 != j))
                                    switch (j) {
                                        case 1: vx = jj*v0x;
                                            vy = jj*v0y;
                                            vz = jj*v0z;
                                            break;
                                        case 3: vx = jj*v1x;
                                            vy = jj*v1y;
                                            vz = jj*v1z;
                                            break;
                                        case 5: vx = jj*v2x;
                                            vy = jj*v2y;
                                            vz = jj*v2z;
                                            break;
                                        case 7: vx = jj*v3x;
                                            vy = jj*v3y;
                                            vz = jj*v3z;
                                            break;
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

    // копируем систему столько раз, сколько потребуется.
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

    // "расширяем" систему до необходимой плотности.
    for (int ii = 0; ii < NP; ii++) {
        particles[ii].x = particles[ii].x*Betta;   // множитель Betta - это коэффициент расширения
        particles[ii].y = particles[ii].y*Betta;   // на него умножаются все координаты и параметры объема
        particles[ii].z = particles[ii].z*Betta;
    }

	double Mx = 0.0;
	for (int i = 0; i < NP; i++) {
		Mx += (particles[i].y + A)*particles[i].vz - (particles[i].z + A)*particles[i].vy;
	}

	printf("\n Mx = %.15le\n", Mx);

    // сохранить систему по старому алгоритму
    // загрузить систему по новому алгоритму
    FILE *save_file = fopen("new.txt", "w+");
    fprintf(save_file, "%d\n", NP);
    fprintf(save_file, "%.15le\n", A);
    fprintf(save_file, "%.15le\n", L);
    for (int i = 0; i < NP; i++) {
        fprintf(save_file, "%.15le %.15le %.15le\n", particles[i].x, particles[i].y, particles[i].z);
        fprintf(save_file, "%.15le %.15le %.15le\n", particles[i].vx, particles[i].vy, particles[i].vz);
    }
    fclose(save_file);

    load_seed("new.txt");
}

/*
 функция изменения состояния частиц в соответствии с наступившим
 в системе событием.
 параметры:
 im - номер частицы
 jm - номер второй частицы или номер границы
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
            int k = im - NP;
            p1.vy = particles[k].vy;
            p1.vz = particles[k].vz;

            if (p1.vx * particles[k].vx < 0) {
                p1.vx = -particles[k].vx;
            }
            else {
                p1.vx = particles[k].vx;
            }
        }
        if (jm >= NP) {
            int k = jm - NP;
            p2.vy = particles[k].vy;
            p2.vz = particles[k].vz;

            if (p2.vx * particles[k].vx < 0) {
                p2.vx = -particles[k].vx;
            }
            else {
                p2.vx = particles[k].vx;
            }
        }

		double E1 = p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz + p2.vx*p2.vx + p2.vy*p2.vy + p2.vz*p2.vz;
		double r = dx*dx + dy*dy + dz*dz;
        /////////////////////////////
		/*
        for (int t = 0; t < particles_for_check_count; t++) {
            if (im == particles_for_check[t] || jm == particles_for_check[t]) {
                FILE *save_file3 = fopen("history.txt", "a");
                fprintf(save_file3, "\n BEFORE collision:");
                fprintf(save_file3, "\n %d vx = %.15le, vy = %.15le, vz = %.15le", im, p1.vx, p1.vy, p1.vz);
                fprintf(save_file3, "\n %d vx = %.15le, vy = %.15le, vz = %.15le\n", jm, p2.vx, p2.vy, p2.vz);
                fclose(save_file3);
            }
        }
		*/
        /////////////////////////////

		//r = dx*dx + dy*dy + dz*dz;
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

		double E2 = p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz + p2.vx*p2.vx + p2.vy*p2.vy + p2.vz*p2.vz;

		if (abs(E1 - E2) > 1.0e-14) {
			if (r == 4.0) {
				printf("\n r = %.15le \n", r);
				printf("\n E1 = %.15le, E2 = %.15le \n", E1, E2);
				printf("\n p1.vx = %.15le, p1.vy = %.15le, p1.vz = %.15le\n", p1.vx, p1.vy, p1.vz);
				printf("\n p2.vx = %.15le, p2.vy = %.15le, p2.vz = %.15le\n", p2.vx, p2.vy, p2.vz);
				printf("\n p1 = %d, p2 = %d, dt = %.15le\n", im, jm, dt);
				printf("a");
			}
		}

        /////////////////////////////
		/*
        for (int t = 0; t < particles_for_check_count; t++) {
            if (im == particles_for_check[t] || jm == particles_for_check[t]) {
                FILE *save_file4 = fopen("history.txt", "a");
                fprintf(save_file4, "AFTER collision:");
                fprintf(save_file4, "\n %d vx = %.15le, vy = %.15le, vz = %.15le", im, p1.vx, p1.vy, p1.vz);
                fprintf(save_file4, "\n %d vx = %.15le, vy = %.15le, vz = %.15le\n", jm, p2.vx, p2.vy, p2.vz);
                fclose(save_file4);
            }
        }
		*/
        /////////////////////////////

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
			int k = im;
			if (particles[im].i_copy < NP)
				k = k - NP;
			destroy_virt_particle(k);
			p1.i_copy = -1;
		}
		if (particles[jm].i_copy > 0) {
			int k = jm;
			if (particles[jm].i_copy < NP)
				k = k - NP;
			destroy_virt_particle(k);
			p2.i_copy = -1;
		}

        create_virt_particle(e1);
        create_virt_particle(e2);
    } else
    if (jm == -1) {
        p1.vx = -p1.vx;
        particles[im] = p1;
    } else {
        if (jm != -100) {
            short end = boxes_yz[p1.y_box][p1.z_box][p1.x_box].end;
            Box p1_box = boxes_yz[p1.y_box][p1.z_box][p1.x_box];

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
                create_virt_particle(im);
                p1 = particles[im];
            }
            p1.box_i = ++boxes_yz[p1.y_box][p1.z_box][p1.x_box].end;
            boxes_yz[p1.y_box][p1.z_box][p1.x_box].particles[p1.box_i] = im;
            particles[im] = p1;
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
 функция "шаг", основной цикл программы
 */
void step() {
	particle p1;
	int i, im, jm;
	double time = 0.0;
	bool need_virt_particle_retime;

	COLL_COUNT = 0;

	while (COLL_COUNT < NP / 2 || jm != -100) {
		im = time_queue[1].im;
		jm = time_queue[1].jm;

		//if (im == 3058 || im == 2826 || jm == 3058 || jm == 2826)
        //    printf("\n Next event: %d %d %le\n", time_queue[1].im, time_queue[1].jm, particles[im].dt);
		iim = im;
		ijm = jm;
        /////////////////////////
		/*
        for (int t = 0; t < particles_for_check_count; t++) {
            if (im == particles_for_check[t] || jm == particles_for_check[t] || particles[im].i_copy == particles_for_check[t] || particles[jm].i_copy == particles_for_check[t]) {
                FILE *save_file = fopen("history.txt", "a");
                fprintf(save_file, "\n Event #%d: %d %d %.15le \n", COLL_COUNT, time_queue[particles[im].ti].im, time_queue[particles[im].ti].jm, particles[im].dt);
                fprintf(save_file, "   %d particle: x %.5le | y %.5le | z %.5le\n", im, particles[im].x, particles[im].y, particles[im].z);
                fprintf(save_file, "   %d particle: vx %.5le | vy %.5le | vz %.5le\n", im, particles[im].vx, particles[im].vy, particles[im].vz);
                fprintf(save_file, " %d i_copy \n", particles[im].i_copy);
                if (particles[im].i_copy >= 0) {
                    fprintf(save_file, "   %d particle: x %.5le | y %.5le | z %.5le\n", particles[im].i_copy, particles[particles[im].i_copy].x, particles[particles[im].i_copy].y, particles[particles[im].i_copy].z);
                    fprintf(save_file, "   %d particle: vx %.5le | vy %.5le | vz %.5le\n", particles[im].i_copy, particles[particles[im].i_copy].vx, particles[particles[im].i_copy].vy, particles[particles[im].i_copy].vz);
                }

                if (jm >= 0) {
                    fprintf(save_file, "   %d particle: x %.5le | y %.5le | z %.5le\n", jm, particles[jm].x, particles[jm].y, particles[jm].z);
                    fprintf(save_file, "   %d particle: vx %.5le | vy %.5le | vz %.5le\n", jm, particles[jm].vx, particles[jm].vy, particles[jm].vz);
                    fprintf(save_file, " %d i_copy \n", particles[jm].i_copy);
                }
                fprintf(save_file, "\n");

                fclose(save_file);
            }
        }
		*/
        /////////////////////////

        delete_event(1);
        particles[im].ti = -1;
        if (jm >= 0) particles[jm].ti = -1;

        p1 = particles[im];

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

        /////////////////////////
		/*
        for (int t = 0; t < particles_for_check_count; t++) {
            if (im == particles_for_check[t] || jm == particles_for_check[t] || particles[im].i_copy == particles_for_check[t] || particles[jm].i_copy == particles_for_check[t]) {
                FILE *save_file = fopen("history.txt", "a");
                fprintf(save_file, "\n After the reform: \n");
                fprintf(save_file, "   %d particle: x %.5le | y %.5le | z %.5le\n", im, particles[im].x, particles[im].y, particles[im].z);
                fprintf(save_file, "   %d particle: vx %.5le | vy %.5le | vz %.5le\n", im, particles[im].vx, particles[im].vy, particles[im].vz);
                fprintf(save_file, " %d i_copy \n", particles[im].i_copy);

                if (jm >= 0) {
                    fprintf(save_file, "\n After the reform: \n");
                    fprintf(save_file, "   %d particle: x %.5le | y %.5le | z %.5le\n", jm, particles[jm].x, particles[jm].y, particles[jm].z);
                    fprintf(save_file, "   %d particle: vx %.5le | vy %.5le | vz %.5le\n", jm, particles[jm].vx, particles[jm].vy, particles[jm].vz);
                    fprintf(save_file, " %d i_copy \n", particles[jm].i_copy);
                }

                fprintf(save_file, "\n");

                fclose(save_file);
            }
        }
		*/
        /////////////////////////
    }

    particle *p = particles;
    for (i = 0; i < NP*2; ++i, ++p)
        (*p).t -= time;
    Event *t = time_queue;
    ++t;
    for (i = 1; i < last; ++i, ++t)
        (*t).t -= time;
    time = 0.0;
}

/*
   Функция получения профиля плотности системы
   Перед снятием характеристик производится синхронизация всех частиц
   Аргументы:
   file_name - имя файла для сохранения данных
 */
void image(int steps, short accuracy, std::string file_name) {
    const int W = (int(L) + 1)*2 * accuracy; // number of dots for all system
    int img[10000];

    printf("INFO: Image started for %d steps with accuracy %d\n", steps, accuracy);
    for (short g = 0; g < 10000; ++g) img[g] = 0;
    for (short h = 0; h < steps; ++h) {
        step();
        // очистка экрана и вывод информации о прогрессе
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("%d / %d", h + 1, steps);

        for (int i = 0; i < NP; ++i) {
            int k = int((L + particles[i].x - particles[i].vx * particles[i].t)*10.0 * accuracy) / 10;
            ++img[k];
        }
    }

    FILE *save_file = fopen(file_name.c_str(), "w+");
    for (int f = 0; f < W; ++f) {
        fprintf(save_file, "%d\n", img[f]);
    }
    fclose(save_file);

    printf("\nINFO: Image completed. Information saved to file:%s\n", file_name.c_str());
}

/*
   Функция создания разреза системы в пике
   позволяет просматривать расположение частиц в слое,
   параллельном идеальной стенке.
   Аргменты:
   x1 - начальная X координата "среза"
   x2 - конечная X координата "среза"
   file_name - имя файла для сохранения данных
 */
void profile(double x1, double x2, int steps, std::string file_name) {
    double x, y, z;
    int i, j;
    FILE *profile_file = fopen(file_name.c_str(), "w+");

    printf("INFO: Profile started.\n");

    for (i = 0; i < steps; ++i) {
        step();
        for (j = 0; j < NP; ++j) {
            x = L + particles[j].x - particles[j].vx * particles[j].t;
            if (x <= x2 && x >= x1) {
                y = A + particles[j].y - particles[j].vy * particles[j].t;
                z = A + particles[j].z - particles[j].vz * particles[j].t;
                // сохраняем Y и Z координаты частицы
                fprintf(profile_file, "%.15le\n", y);
                fprintf(profile_file, "%.15le\n", z);
            }
        }
    }

    fclose(profile_file);

    printf("INFO: Profile completed. Information saved to file: %s\n", file_name.c_str());
}

void compress(double compress_to_etta) {
	printf("INFO: Start to change system density...\n");
	double etta = (PI * NP) / (6.0 * A * A * (L - 1.0));
	double min = 1.0e+100, max = -1.0, x;
	// задаём плотность с точностью в 12 знаков
	while (fabs(etta - compress_to_etta) > 1.0e-12) {
		if (etta < compress_to_etta) {
			max = -100.0;
			min = 1.1e+10;

			for (int i = 0; i < NP; ++i) {
				x = L + particles[i].x - particles[i].vx * particles[i].t;
				if (x < min) min = x;
				if (x > max) max = x;
			}

			min = min - 1.0;
			max = 2.0 * L - max - 1.0;
			if (max < min) min = max;

			// сжимаем не впритык к частицам и не слишком быстро
			min = min / 1.1;
			if (min < 0.1e-8) min = 0.01e-30;
			if (min > 0.01) min = 0.01;
		}
		else {
			double L_ideal = ((PI * NP) / compress_to_etta) / (6.0 * A * A) + 1;
			min = L - L_ideal;
			if (min < -0.01) min = -0.01; // шаг расширения системы
		}
		// изменяем систему
		L -= min;

		save("tmp");
		load_seed("tmp");

		for (short i = 0; i < 10; ++i)
			step();
		etta = (PI * NP) / (6.0 * A * A * (L - 1.0));

		// очистка экрана и вывод информации о прогрессе
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("etta = %.15le, should be equal to %.15le", etta, compress_to_etta);
	}
	printf("\n INFO: System density was sucessfully changed to %.15le \n", etta);
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
            int NN = 8;
            double etta = 0.3;
			K = 12;
			K2 = 32;

            command_file >> NN;
            command_file >> etta;
            command_file.getline(parameter, 255, '\n');
            new_seed(NN, etta);
            etta = (PI * NP) / (6.0 * A * A * (L - 1.0));
			printf("| A    | %.15le |\n| L    | %.15le |\n| N    | %d                  |\n| etta | %.15le |\n", 2.0 * A, 2.0 * L, NP, etta);
        }
        if (str_command.compare("load") == 0) {
            command_file.getline(parameter, 255, '\n');
            load_seed(parameter);
            double etta = (PI * NP) / (6.0 * A * A * (L - 1.0));
			printf(" A = %.15le; \n L = %.15le; \n etta = %.15le; \n", 2.0 * A, 2.0 * L, etta);
        }
        if (str_command.compare("step") == 0) {
            command_file >> steps;
            command_file.getline(parameter, 255, '\n');

            printf(" INFO: Step start. \n");

            start = clock();

            for (i = 0; i < steps; ++i) {
                step();
				check_particles();

                printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                printf("%d / %d", i + 1, steps);
            }

            end = clock();
            result = end - start;

            printf("\n INFO: finished %d collisions per particle \n", steps);
            printf("Total Time = %f seconds. \n", double(result / CLOCKS_PER_SEC));
        }
        if (str_command.compare("image") == 0) {
            command_file >> steps; // число соударений на одну частицу за время измерений
            command_file >> i; // точность. число точек графика на один радиус частицы
            command_file.getline(parameter, 255, '\n');
            image(steps, i, parameter);
        }
        if (str_command.compare("profile") == 0) {
            double x1, x2;
            command_file >> x1;
            command_file >> x2;
            command_file >> steps;
            command_file.getline(parameter, 255, '\n');
            profile(x1, x2, steps, parameter);
        }
        if (str_command.compare("save") == 0) {
            command_file.getline(parameter, 255, '\n');
            save(parameter);
			printf("\n INFO: particles coordinates saved to '%s' \n", parameter);
        }
        if (str_command.compare("compress") == 0) {
            double etta;
            command_file >> etta; // требуемая плотность
            command_file.getline(parameter, 255, '\n');
            compress(etta);
        }
        if (str_command.empty()) break;
    }
    FILE *result_flag = fopen("result", "w+");
    fclose(result_flag);
}


int main(array<System::String ^> ^args)
{
    //K = 3;
    //K2 = 8;
    //new_seed(2, 0.3);

	//K = 12;
	//K2 = 32;
	//new_seed(8, 0.3);

    FILE *save_file = fopen("history.txt", "w");
    fclose(save_file);
    //particles_for_check_count = 0;
    //particles_for_check[0] = 5154;
	//particles_for_check[1] = 3165;

	init("program.txt");
    return 0;
}