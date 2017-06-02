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
ИМПОРТИРУЕМЫЕ МОДУЛИ
*/
#include "stdafx.h"  // стандартная библиотека
#include <time.h>    // библиотека для подстчёта времени
#include <fstream>   // библиотека для работы с файлами
#include <cmath>     // библиотека для вызова математических функций


/*
ОБЪЯВЛЕНИЕ ПЕРЕМЕНЫХ
*/

// Количество ячеек по Y, Z и X(K2), на которые разбивается объём.
// Количество ячеек расчитывается динамически при загрузке системы
// или при посеве (см. функцию load).
long int K, K2;

// число частиц во всём объёме, переопределяется в функциях load и new_seed
long int NP;

// глобальный счётчик столкновений в системе
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

// параметры объёма, задаются в load()
double A, dA, L, dL;

// Глобальная переменная для подсчёта общей кинетической энергии всех частиц
double global_E = 0.0;

double dV_left_wall, dV_right_wall, N_collissions_left_wall, N_collissions_right_wall, dP_time;

// индекс последнего элемента в очереди событий.
long int last;

// объект "событие"
typedef struct Event_ {
	double t;
	long int im, jm, vp1, vp2;
} Event;

// очередь событий - оптимально 16384 элемента
// (это должно быть число-степень двойки, большее чем максимальное число частиц)
//Event *time_queue = new Event[262144];
Event *time_queue = new Event[16384];

// объект "частица"
// x, y, z - координаты частицы
// vx, vy, vz - проекции скоростей частицы
// t - собственное время частицы
// dt - время до ближайшего события этой частицы
// max_x - максимальная координата Х до которой долетает
//         частица при диффузном отражении
// x_box, y_box, z_box - номер ячейки, в которой находится частица
// ti - номер события частицы в дереве событий
// box_i - номер частицы в ячейке
typedef struct particle_ {
	double x[4], y[4], z[4], vx, vy, vz, t, dt, x_max;
	long int x_box[4], y_box[4], z_box[4], ti, box_i[4];
} particle;

// массив частиц, размер массива N + округление
// в большую сторону к числу дающее степень двойки.
//particle *particles = new particle[262144];
particle *particles = new particle[16384];

// клетка. Объём системы разделён на множество клеток,
// каждая клетка содержит в себе несколько виртуальных частиц
// x1, y1, z1, x2, y2, z2 - координаты конца и начала каждой ячейки
// particles[] - список всех частиц, находящихся в данной ячейке
// end - индекс последней частицы в списке частиц данной ячейки
typedef struct Box_ {
	double x1, y1, z1, x2, y2, z2;
	long int particles[32];
	int end;
} Box;

// массив клеток для всего объёма
Box boxes_yz[100][100][1000];


// массив с номерами частиц, для которых надо сохранять историю событий
// используется на случай отладки программы для сохранения истории событий
// выбранных частиц
long int particles_for_check[100];
long int particles_for_check_count = 0;

double *angles = new double[200];

/*
Эта функция выводит на экран параметры текущей системы:
A - размер системы (области объёма, где может находиться центр частицы) по Y и Z.
L
N - число частиц в системе
etta - средняя относительная плотность системы
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
Функция возвращает наибольшее собственное время частиц в системе,
что позволяет синхронизовать все частицы по времени
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
Функция для проверки состояния системы, в ней мы проверяем, что
все частицы находятся внутри системы, для каждой частицы у нас рассчитано
ближайшее событие, частицы находятся в правильных ячейках системы и пр.

В случае возникновения любых проблем данная функция прекращает работу программы
и выводит дополнительную информацию об обнаруженной проблеме.

Функция работает медленно, необходимо использовать в целях проверки
изменений в программе.
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

		// проверяем индексы ячеек для всех частиц
		if ((p1.x_box[0] > K2) || (p1.y_box[0] > K) || (p1.z_box[0] > K) ||
			(p1.x_box[0] < 0) || (p1.y_box[0] < 0) || (p1.z_box[0] < 0)) {
			throw "Particle locates in incorrect cell.";
		}

		// проверяем что частицы находятся в объёме
		if ((p1.x[0] > L && p1.x[0] < 0.0) ||
			(p1.y[0] > A && p1.y[0] < 0.0) ||
			(p1.z[0] > A && p1.z[0] < 0.0)) {
			printf(" \n Particle %ld, %.15le, %.15le, %.15le \n ", i, p1.x[0], p1.y[0], p1.z[0]);
			throw "Particle is out of the system boundaries.";
		}

		// проверяем что частицы находятся в правильных ячейках
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

		// Проверяем что частица записана в одную из ячеек в системе
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
		// Проверяем на правильное ли событие в линейке времён ссылается частица
		if (time_queue[p1.ti].im != i && time_queue[p1.ti].jm != i) {
		Event e = time_queue[p1.ti];
		printf("\n i = %d ; im = %d ; jm = %d ; ti = %d ", i, e.im, e.jm, p1.ti);
		throw "Particle has no correct link to the event.";
		}
		*/
	}

	// Проверяем текущее значение глобальной кинетической энергии системы со
	// значением энергии, которое было при загрузке системы из файла в load
	if (abs(E - global_E) > global_E*0.1e-8) {
		printf("\nENERGY was changed: \n E_seed = %.15le \n E_now= %.15le \n", global_E, E);
		throw "ENERGY was changed.";
	}

	return 0;
}


/*
Функция подъёма элемента по очереди событий к началу очереди

Аргументы:
i - позиция, на которой находится элемент в данный момент
t - время до наступления данного события
*/
int get_up(long int i, double &t) {
	long int j = i >> 1;  // это сдвиг вправо, то же самое что j = i/2, только быстрее
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
Функция добавления события в очередь событий

Аргументы:
i - частица, с участием которой произойдёт новое событие
j - номер частицы с которой столкнётся частица i или номер события
соударения со стенкой или прохождения через периодические
граничные условия или между ячейками системы
*/
void add_event(long int &i, long int &j, long int &vp1, long int &vp2) {
	/*
	Рассчитываем полное время нового события от начала отсчёта
	глобального времени системы, таким образом мы получаем время t,
	сравнивая которое мы можем определить какое из событий в системе
	произойдёт раньше
	*/
	double t = particles[i].t + particles[i].dt;

	/*
	Находим позицию в дереве времён для нового события.
	Изначально помещаем это событие вниз дерева и позволяем
	ему подняться по дереву, если данное событие произойдёт раньше чем
	другие события
	*/
	particles[i].ti = get_up(last, t);

	/*
	Если новое событие - это событие столкновения двух частиц, то
	необходимо для второй частицы сохранить данные о её новом событии
	*/
	if (j >= 0) {
		particles[j].dt = particles[i].dt;
		particles[j].ti = particles[i].ti;
	}

	// записываем новое событие в выбранную ячейку в дереве времён
	time_queue[particles[i].ti].im = i;
	time_queue[particles[i].ti].jm = j;
	time_queue[particles[i].ti].vp1 = vp1;
	time_queue[particles[i].ti].vp2 = vp2;
	time_queue[particles[i].ti].t = t;

	// увеличиваем число событий на 1
	last++;
}


/*
Функция удаления события из очереди событий

Аргументы:
i - позиция удаляемого элемента в очереди событий
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
Функция удаления событий частицы из линейки событий

Аргументы:
i - номер частицы или образа, события которого необходимо удалить
*/
void clear_particle_events(long int &i) {
	// particles[i].ti - ссылка на индекс события, которое хранит сама частица
	int f = particles[i].ti;

	if (f > 0) {
		long int e = -100, pv1 = 0, pv2 = -1;
		long int kim = time_queue[f].im;
		long int kjm = time_queue[f].jm;

		if (time_queue[f].im == i) {
			/*
			Если мы удаляем событие столкновения двух частиц, то
			необходимо переместить вторую частицу в то же время,
			в котором находится частица, для которой мы удаляем
			это событие и заменить событие столкновения двух частиц
			на событие "-100", когда вторая частица просто долетит до
			места предполагаемого столкновения и после этого для неё
			будет рассчитано новое событие.
			*/
			if (time_queue[f].jm >= 0) {
				double dt = particles[kim].t - particles[kjm].t;  // разница в текущем времени

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
			double dt = particles[kjm].t - particles[kim].t;  // разница в текущем времени

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
Функция расчёта ближайшего события для частицы

Аргументы:
i - номер частицы, для которой мы должны рассчитать ближайшее событие
*/
void retime(long int &i) {
	particle p1 = particles[i];
	Box p1_box = boxes_yz[p1.y_box[0]][p1.z_box[0]][p1.x_box[0]];
	long int jm = -1;  // переменная для сохранения типа ближайшего события
	double dt, dt_min = 1.0e+20;  // переменные для рассчёта времени для ближайшего события

	if (p1.vx < 0.0) {
		dt_min = (p1_box.x1 - p1.x[0]) / p1.vx;
		jm = -2;  // событие пересечения границы Х1 ячейки, в которой находится частица 

				  // если мы находимся вблизи идеальной стенки, то рассчитать время соударения
				  // с идеальной стенкой
		if (p1.x_box[0] == 1) {
			dt_min = (p1_box.x1 + particle_R - p1.x[0]) / p1.vx;
			jm = -1;  // событие столкновения с идеальной стенкой
		}
	}
	else {
		dt_min = (p1_box.x2 - p1.x[0]) / p1.vx;
		jm = -4;  // событие пересечения границы Х2 ячейки, в которой находится частица

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
			// если мы находимся вблизи идеальной стенки, то рассчитать время соударения
			// с идеальной стенкой
			if (p1.x_box[0] == K2 - 1) {
				dt_min = (p1_box.x2 - particle_R - p1.x[0]) / p1.vx;
				jm = -1;  // событие столкновения с идеальной стенкой
			}
		}
	}

	if (p1.vy < 0.0) {
		dt = (p1_box.y1 - p1.y[0]) / p1.vy;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -5;  // событие пересечения границы Y1 ячейки, в которой находится частица  
		}
		if (p1.y_box[0] == 1) {
			dt = (p1_box.y1 + particle_Rx2 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -12;  // событие рождения образа частицы
			}
		}
		if (p1.y_box[0] == K - 1) {
			dt = (p1_box.y2 - particle_Rx2 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -11;  // событие удаления образа частицы
			}
		}
	}
	else {
		dt = (p1_box.y2 - p1.y[0]) / p1.vy;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -6;  // событие пересечения границы Y2 ячейки, в которой находится частица
		}
		if (p1.y_box[0] == K - 1) {
			dt = (p1_box.y2 - particle_Rx2 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -11;  // событие рождения образа частицы
			}
		}
		if (p1.y_box[0] == 1) {
			dt = (p1_box.y1 + particle_Rx2 - p1.y[0]) / p1.vy;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -12;  // событие удаления образа частицы
			}
		}
	}

	if (p1.vz < 0.0) {
		dt = (p1_box.z1 - p1.z[0]) / p1.vz;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -7;  // событие пересечения границы Z1 ячейки, в которой находится частица
		}
		if (p1.z_box[0] == 1) {
			dt = (p1_box.z1 + particle_Rx2 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -14;  // событие рождения образа частицы
			}
		}
		if (p1.z_box[0] == K - 1) {
			dt = (p1_box.z2 - particle_Rx2 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -13;  // событие удаления образа частицы
			}
		}
	}
	else {
		dt = (p1_box.z2 - p1.z[0]) / p1.vz;
		if (dt < dt_min) {
			dt_min = dt;
			jm = -8;  // событие пересечения границы Z2 ячейки, в которой находится частица
		}
		if (p1.z_box[0] == K - 1) {
			dt = (p1_box.z2 - particle_Rx2 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -13;  // событие рождения образа частицы
			}
		}
		if (p1.z_box[0] == 1) {
			dt = (p1_box.z1 + particle_Rx2 - p1.z[0]) / p1.vz;
			if ((dt > 0) && (dt < dt_min)) {
				dt_min = dt;
				jm = -14;  // событие удаления образа частицы
			}
		}
	}

	double temp, dx, dy, dz, dvx, dvy, dvz, d, dv, bij, dr;
	long int s, n, r, q, w, vp1 = 0, vp2 = -1;

	for (long int iv = 0; iv < 4; iv++) {
		if (p1.x[iv] > 0.0) {
			// Проходим по ячейкам, ближайшим к ячейке, в которой находится частица i
			for (r = p1.x_box[iv] - 1; r < p1.x_box[iv] + 2; r++)
				for (q = p1.y_box[iv] - 1; q < p1.y_box[iv] + 2; q++)
					for (w = p1.z_box[iv] - 1; w < p1.z_box[iv] + 2; w++) {

						// если индекс ячейки выходит за границу системы то
						// переходим на следующий шаг цикла
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

						// Проходим по всем частицам в выбранной ячейке и проверяем возможность
						// столкновения этих частиц с частицей i
						for (s = 0; s <= boxes_yz[q][w][r].end; ++s) {
							n = boxes_yz[q][w][r].particles[s];
							long int jv = 0;

							if (n >= NP) {
								jv = n / NP; // номер изначальной частицы
								n = n % NP;  // номер образа
							}

							/*
							Cохраняем в переменную p все данные частицы n,
							чтобы далее записывать все операции короче и без обращения
							в глобальную память
							*/
							particle p = particles[n];

							/*
							Рассчитываем разницу собственного времени двух частиц,
							это необходимо чтобы синхронизовать их между собой и рассчитывать
							возможность и время столкновения в системе, где обе частицы будут
							иметь одинаковое собственное время
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
								// рассчитываем дискриминант в уравнении для вычисления времени соударения
								d = bij * bij + dv * dr;

								// если дискриминант больше нуля то соударение возможно
								if (d > 0.0) {
									dt = -(sqrt(d) + bij) / dv;

									/*
									К разнице в собственном времени частиц прибавляем время до их соударения,
									в результате получаем время dt через которое это соударение произойдёт для
									частицы p, таким образом мы получаем возможность сравнить время только что
									рассчитанного соударения со временем ближайшего события для частицы p и
									принять решение - должны ли мы перезаписать событие для частицы p или нет.
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
	Если новое событие - это столкновение двух частиц, то для второй частицы
	необходимо перезаписать её предыдущее ближайшее событие и синхронизовать
	эту частицу по времени с частицей i
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
	В случае если при расчёте ближайшего события мы получили отрицательное время
	необходимо прервать выполнение программы и распечатать отладочную информацию
	о рассчитанном событии
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
Функция обмена частицы и "образа" при пересечении частицей периодических
граничных условий.

Аргументы:
im - номер частицы, пересекающей периодические граничные условия
jm - номер границы, через которую проходит центр частицы
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
Функция создания нового "образа".

Аргументы:
i - номер частицы, образ которой необходимо создать
need_to_check - флаг, позволяющий не проверять нужен ли образ для данной
частицы или нет, и сразу создавать образ (когда мы уверены,
что образ необходимо создать).
*/
void create_virt_particle(long int &im, long int &jm) {
	Box box1;

	switch (jm)
	{
	case -11:     // частица на границе Y = А-2R
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
	case -12:     // частица на границе Y = 2R
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
	case -13:     // частица на границе Z = A-2R
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
	case -14:     // частица на границе Z = 2R
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

	// если у частицы есть два образа, то необходимо создать третий образ
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
Функция загрузки информации о некоторой частице из файла.
При загрузке системы мы последовательно считываем информацию
о каждой частице с помощью этой функции.

Аргументы:
loading_file - ссылка на файл, из которого необходимо считать
данные о частице или образе.
*/
void load_information_about_one_particle(FILE *loading_file) {
	double a1, a2, a3;
	long i, x_box, y_box, z_box, end;

	// Считывем номер частицы
	fscanf(loading_file, "%ld", &i);

	// Считываем координаты частицы
	fscanf(loading_file, "%le %le %le", &a1, &a2, &a3);
	particles[i].x[0] = a1;
	particles[i].y[0] = a2;
	particles[i].z[0] = a3;

	// Считываем информацию о скоростях частицы
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
			// Определяем в какую ячейку необходимо записать новую частицу
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

			// записываем частицу в ячейку
			end = particles[i].box_i[j] = ++boxes_yz[y_box][z_box][x_box].end;
			boxes_yz[y_box][z_box][x_box].particles[end] = i + NP*j;

			/*
			Если в одной ячейке находится слишком много частиц,
			то завершить программу, выведя сообщение об ошибке.
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
	Проверяем что в данной ячейке нет повторяющихся номеров частиц
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
Функция загрузки данных о системе из файла.
Загружаются параметры объёма, координаты и скорости всех частиц.

Аргументы:
file_name - имя файла, содержащего все данные о системе.
*/
void load_seed(std::string file_name) {
	double a1, x, y, z;
	long h;
	FILE *loading_file;

	global_E = 0.0;

	/*
	Открываем указанный файл для чтения,
	пишем на экран информацию об ошибке и выходим из программы если
	указанного файла не существует или его невозможно открыть для чтения.
	*/
	try {
		loading_file = fopen(file_name.c_str(), "r");
	}
	catch (...) {
		printf("\n ERROR: Can't read file '%s' - is it exist? \n", file_name.c_str());
		exit(1);
	}

	// считываем плотность
	fscanf(loading_file, "%le", &a1);

	// считываем количество частиц в системе
	fscanf(loading_file, "%ld", &h);
	NP = h;
	// считываем значение A
	fscanf(loading_file, "%le", &a1);
	A = a1;
	// считываем значение L 
	fscanf(loading_file, "%le", &a1);
	L = a1;
	// считываем значение R 
	//fscanf(loading_file, "%le", &a1);
	//change_R(a1);

	// EN: search for the appropriate values for dA, dL, K, K2
	// RU: ищем подходящее значение для dA, dL, K, K2
	dA = 2.2*particle_R;
	dL = 2.2*particle_R;
	K = long(A / dA) + 1;  // количество ячеек по Y и Z
	K2 = long(L / dL) + 1; // количество ячеек по X
	dA = A / (K - 1);       // точные размеры ячеек по X, Y и Z
	dL = L / (K2 - 1);      // выбранные так, чтобы, ячейки совпадали
							// с размерами объёма

							/*
							После того как мы рассчитали количество ячеек и их размеры
							необходимо инициализировать данные каждой ячейки, указать
							границы каждой ячейки по X, Y, Z, т.к. эти данные будут
							использоваться для рассчёта времени до пересечения границ
							ячеек частицами, которые находятся в этих ячейках.
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
	Обнуляем все данные для всех частиц и образов перед считыванием
	файла с данными.
	Для образов, которые не описаны в файле сохранения, останутся эти
	начальные значения.
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
	Загружаем данные о всех частицах и "образах" из файла сохранения
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
	Инициализируем очередь событий пустыми событиями и устанавливаем
	указатель конца списка на его начало.
	*/
	last = 1;
	time_queue[0].t = 0.0;
	time_queue[0].im = -1;
	for (long i = 1; i < NP; i++) time_queue[i].t = 1.0E+20;

	/*
	Рассчитываем новые события для всех частиц и существующих в системе "образов"
	*/
	for (long i = 0; i < NP; i++) {
		retime(i);
	}
}


/*
Функция сохранения состояния системы в текстовый файл.
перед сохранением производится синхронизация частиц по времени.

Аргументы:
file_name - имя файла, в который будет записана информация
*/
void save(std::string file_name) {
	double x, y, z, dt;
	double t_global = get_maximum_particle_time();

	FILE *save_file = fopen(file_name.c_str(), "w+");

	double etta = (4.0 * PI * particle_R3 * NP) / (3.0 * A * A * (L - particle_Rx2));
	fprintf(save_file, "%.15le\n", etta);

	// RU: сохраняем информацию о количестве частиц и размерах системы
	// EN: save information about count of particles and size of the system
	fprintf(save_file, "%ld\n", NP);
	fprintf(save_file, "%.15le\n", A);
	fprintf(save_file, "%.15le\n", L);
	//fprintf(save_file, "%.15le\n", particle_R);

	// RU: сохраняем координаты всех частиц и их скорости: x, y, z, vx, vy, vz
	// EN: we need to save all coordinates of particles: x, y, z, vx, vy, vz
	for (long i = 0; i < NP; ++i) {
		fprintf(save_file, "%ld\n", i);

		/*
		Синхронизируем частицу i с глобальным временем системы,
		в итоге данные всех частиц будут записаны для одного момента времени,
		который совпадает со временем последнего произошедшего в системе события.
		*/
		dt = t_global - particles[i].t;

		/*
		При сохранении состояния системы мы указываем
		абсолютные координаты частиц и параметры системы, так, чтобы
		центр системы был в точке (L/2, A/2), а не в (0, 0).
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
Процедура нового "посева" системы.

Аргументы:
NN - число частиц в ребре объёмо центрированного кристалла, на основе
которого делается изначальный посев
etta - начальная плотность, которую необходимо задать в системе
lamda - расстояние между узлами
*/
void new_seed(long particles_count, double etta) {
	long i = 0, j, j1, j2, j3;

	// default 16x48
	double lamda = 2.0*particle_R + 0.1*particle_R;
	long Qyz = particles_count;  // число узлов по у и z
	long Qx = particles_count * 4; // число узлов по x
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

	//srand(time(NULL));  - совсем случайно
	srand(10);   // - не случайно

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
	Cохраняем и загружаем систему из файла чтобы инициализировать
	все необходимые переменные.
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
			Стираем частицу из ячейки в которой она находилась,
			вместо неё вставляем частицу из конца списка частиц в ячейке и
			уменьшаем число частиц в ячейке на 1.
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
Функция изменения состояния частиц в соответствии с наступившим
в системе событием.

Аргументы:
im - номер частицы
jm - номер второй частицы или номер границы
*/
void reform(long &im, long &jm, long &vp1, long &vp2) {
	particle p1 = particles[im];
	double q1, q2, z, dx, dy, dz;

	/*
	Если jm >= 0 то наступившее событие - это соударение двух
	частиц или образов.
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
		case -1:   // соударение с идеальной стенкой
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
		case -11:  // пересечение границы Y = A - 2R
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
		case -12:  // пересечение границы Y = 2R
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
		case -13:    // пересечение границы Z = A - 2R
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
		case -14:  // пересечение границы Z = 2R
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
					Стираем частицу из ячейки в которой она находилась,
					вместо неё вставляем частицу из конца списка частиц в ячейке и
					уменьшаем число частиц в ячейке на 1.
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
					После того, как индекс ячейки для частицы был изменён,
					добавляем частицу в новую ячейку (или возвращаем её в старую
					ячейку, если индекс ячейки не изменился).
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
Функция "шаг", основной цикл программы,
производит рассчёт динамики системы в течении 1 соударения
на каждую частицу  - (NP/2) столкновений в системе.
*/
void step() {
	particle p1;
	long int i, im, jm, vp1, vp2;
	double time = get_maximum_particle_time();

	COLL_COUNT = 0;
	jm = 0;

	while (COLL_COUNT < NP / 2 || jm < 0) {
		// считываем первое событие из линейки событий
		im = time_queue[1].im;
		jm = time_queue[1].jm;
		vp1 = time_queue[1].vp1;
		vp2 = time_queue[1].vp2;

		// удаляем первое событие
		delete_event(1);

		p1 = particles[im];

		/*
		Обнуляем указатель на событие частиц, учавствующих в этом
		событии, приравнивая их -1, так, чтобы они не указывали на
		существующие события.
		*/
		p1.ti = -1;
		if (jm >= 0) particles[jm].ti = -1;

		/*
		Если наступило событие -100, значит просто
		синхронизируем частицу с которой произошло
		данное событие с глобальным временем системы.
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

		// Производим изменения в системе согласно произошедшему событию
		reform(im, jm, vp1, vp2);

		retime(im);

		/*
		Если jm > 0, то это значит, что произошло
		соударение частиц в системе, увеличиваем
		счётчик соударений
		*/
		if (jm >= 0) {
			++COLL_COUNT;
			retime(jm);
		}
	}

	/*
	Запускаем глобальную синхронизацию частиц по собственному времени
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
		// считываем первое событие из линейки событий
		im = time_queue[1].im;
		jm = time_queue[1].jm;
		vp1 = time_queue[1].vp1;
		vp2 = time_queue[1].vp2;

		// удаляем первое событие
		delete_event(1);

		p1 = particles[im];

		/*
		Обнуляем указатель на событие частиц, учавствующих в этом
		событии, приравнивая их -1, так, чтобы они не указывали на
		существующие события.
		*/
		p1.ti = -1;
		if (jm >= 0) particles[jm].ti = -1;

		/*
		Если наступило событие -100, значит просто
		синхронизируем частицу с которой произошло
		данное событие с глобальным временем системы.
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

		// Производим изменения в системе согласно произошедшему событию
		reform(im, jm, vp1, vp2);

		retime(im);

		/*
		Если jm > 0, то это значит, что произошло
		соударение частиц в системе, увеличиваем
		счётчик соударений
		*/
		if (jm >= 0) {
			++COLL_COUNT;
			retime(jm);
		}
	}

	/*
	Запускаем глобальную синхронизацию частиц по собственному времени
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
Функция получения профиля плотности системы

Перед снятием характеристик производится синхронизация всех частиц
с глобальным временем системы.

Аргументы:
steps - число замеров плотности для построения данных о профиле плотности
accuracy - число ячеек по Х на которые разбивается профиль плотности
и в которых считается количество частиц
file_name - имя файла для сохранения данных
*/
void image(long steps, long accuracy, std::string file_name) {
	// число ячеек, в которых будем рассчитывать количество частиц
	const unsigned W = unsigned(L*accuracy);

	printf("\n W = %ld \n", W);

	// массив, в который будем собирать данные по количеству частиц
	long *img = new long[W];
	double t_global, dt;

	printf("INFO: Image started for %ld steps with accuracy %ld small cells per 1.0\n", steps, accuracy);

	// заполняем массив нулями, прежде чем считать количество частиц
	for (unsigned g = 0; g < W; g++) img[g] = 0;

	for (long h = 0; h < steps; h++) {
		step();  // далем одно соударение на частицу между замерами

				 // очистка экрана и вывод информации о прогрессе
		printf("\r%76c\r", ' ');
		printf("%ld / %ld", long(h + 1), steps);

		/*
		Рассчитываем глобальное время системы, которое равно
		самому большому собственному времени частиц в системе.
		*/
		t_global = get_maximum_particle_time();

		for (long i = 0; i < NP; ++i) {
			/*
			Синхронизируем частицу с глобальным временем системы
			*/
			dt = t_global - particles[i].t;

			/*
			Определяем в какую из ячеек по X попадает данная частица
			после синхронизации её по времени с глобальным временем системы
			*/
			long m = long((particles[i].x[0] + particles[i].vx * dt) * 10000.0 * accuracy) / 10000;
			img[m]++;  // увеличиваем количество частиц в ячейке на 1.
		}
	}

	double x = 0.0, delta_x = 1.0 / accuracy;
	FILE *profile_file = fopen(file_name.c_str(), "w+");
	for (unsigned f = 0; f < W; f++) {
		/*
		Сохраняем начальную координату X для ячейки и число частиц,
		которые находились в данном "слое" во время измерений
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

		// очистить строку
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

// очистить строку
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

	// открываем файл на запись
	FILE *profile_file = fopen(file_name.c_str(), "w+");

	printf("INFO: Profile started: %.5le - %.5le.\n", x1, x2);

	double t_global = get_maximum_particle_time();

	particles_count = 0;
	for (j = 0; j < NP; j++) {
		/*
		Синхронизируем частицу по времени с глобальным
		временем системы и проверяем её координату Х,
		если частица находится в указанном диапазоне,
		то сохраняем её Y и Z координаты в файл.
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
		Между измерениями данных о положениях частиц
		делаем указанное количество соударений на частицу в системе
		*/
		for (u = 0; u < steps; u++) {
			step();

			// очистить строку
			printf("\r%76c\r", ' ');
			printf("%ld / %ld", i*steps + u, steps*dots_per_particle);

			t_global = get_maximum_particle_time();

			w = 0;
			for (j = 0; j < particles_count; j++) {
				/*
				Синхронизируем частицу по времени с глобальным
				временем системы и проверяем её координату Х,
				если частица находится в указанном диапазоне,
				то сохраняем её Y и Z координаты в файл.
				*/
				dt = t_global - particles[f[j]].t;
				x = particles[f[j]].x[0] + particles[f[j]].vx * dt;

				y[w] += particles[f[j]].y[0] + particles[f[j]].vy * dt;
				z[w] += particles[f[j]].z[0] + particles[f[j]].vz * dt;

				w++;
			}
		}

		for (m = 0; m < w; m++) {
			// сохраняем Y и Z координаты частицы
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
Функция создания разреза системы в пике,
позволяет просматривать расположение частиц в слое,
параллельном идеальной стенке.

Аргументы:
x1 - начальная X координата "среза"
x2 - конечная X координата "среза"
dots_per_particle - количество снятий данных через равное количество соударений в системе
steps - число соударений между двумя снятиями данных о положении частиц в выбранном слое
file_name - имя файла для сохранения данных
*/
void profile(double x1, double x2, long int dots_per_particle, long int steps, std::string file_name) {
	double x, y, z, dt;
	long int i, j;

	// открываем файл на запись
	FILE *profile_file = fopen(file_name.c_str(), "w+");

	printf("INFO: Profile started: %.5le - %.5le.\n", x1, x2);

	for (i = 0; i < dots_per_particle; ++i) {

		// очистить строку
		printf("\r%76c\r", ' ');
		printf("%ld / %ld", i, dots_per_particle);

		/*
		Между измерениями данных о положениях частиц
		делаем указанное количество соударений на частицу в системе
		*/
		for (j = 0; j < steps; j++)
			step();

		double t_global = get_maximum_particle_time();

		for (j = 0; j < NP; j++) {
			/*
			Синхронизируем частицу по времени с глобальным
			временем системы и проверяем её координату Х,
			если частица находится в указанном диапазоне,
			то сохраняем её Y и Z координаты в файл.
			*/
			dt = t_global - particles[j].t;
			x = particles[j].x[0] + particles[j].vx * dt;

			if (x <= x2 && x >= x1) {
				y = particles[j].y[0] + particles[j].vy * dt;
				z = particles[j].z[0] + particles[j].vz * dt;

				// сохраняем Y и Z координаты частицы
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
x1, x2 - X координаты двух плоскостей, ограничивающих слой кристалла
m1 - количество измерений, которые буду сделаны в течении одной итерации
m2 - количество измерений, которые делаются для нахождения среднего положения центра частицы
m3 - количество итераций измерений
m4 - количество соударений в среднем на одну частицу в системе, которые должны проиойти
между двумя ближайшими измерениями.
*/
void F1(double x1, double x2, int m1, int m2, int m3, int m4, std::string file_name) {

	double t_global, dt;
	particle p1;

	FILE *f1_file = fopen(file_name.c_str(), "w+");
	FILE *debug_file = fopen("debug.log", "w+");

	// Динамически создаём двумерные массивы (это необходимо делать динамически,
	// так как они большого размера и мы не можем создать их в стеке,
	// к тому же, размер массива зависит от значений переменной m1)
	// 200 - по количеству частиц в слое, с запасом (частиц будет 100+, но меньше 200)
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

		// Число частиц, которые будут выбраны для наблюдения в данном слое:
		particles_count = 0;

		// Число частиц, центры которых расположены вдоль одной из осей симметрии слоя:
		particles_count_alfa = 0;
		// Число частиц, центры которых расположены вдоль еще одной оси симметрии слоя
		// данная ось пересекается с первой в центре одной из частиц):
		particles_count_betta = 0;

		// Заполняем все массивы нулями:
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

		// Не рассматриваем частицы, находящиеся вблизи периодических
		// граничных условий, а так же рассматриваем только частицы,
		// находящиеся в данном слое.
		// В этом цикле мы сначала перемещаем все частицы во время системы,
		// а затем проверяем удовлетворяет ли текущее положение центра частицы критериям,
		// чтобы найти все частицы, положение которых необходимо отслеживать:
		t_global = get_maximum_particle_time();
		for (long i = 0; i < NP; i++) {
			p1 = particles[i];
			dt = t_global - p1.t;
			p1.x[0] = p1.x[0] + p1.vx * dt;
			p1.y[0] = p1.y[0] + p1.vy * dt;
			p1.z[0] = p1.z[0] + p1.vz * dt;

			// Отбираем только частицы в данном слое такие, которые находятся на расстоянии
			// больше 6 радиусов от периодических граничных условий:
			if (p1.x[0] >= x1 && p1.x[0] <= x2 &&
				p1.y[0] >= 6.0 * particle_R && p1.y[0] <= A - 6.0 * particle_R &&
				p1.z[0] >= 6.0 * particle_R && p1.z[0] <= A - 6.0 * particle_R) {
				list_of_particles[particles_count] = i;
				particles_count++;
			}
		}

		for (long i = 0; i < m1; i++) {
			// Выводим на экран текущий прогресс расчета:
			printf("\r%76c\r", ' ');
			printf("%ld / %ld", iter*m1 + i + 1, m1*m3);

			// Делаем несколько соударений на частицу между измерениями
			// чтобы положения частиц были достаточно случайны.
			for (int r = 0; r < m4; r++) {
				step();
			}

			t_global = get_maximum_particle_time();
			for (int m = 0; m < particles_count; m++) {
				w = list_of_particles[m];

				// Cинхронизируем время частицы с временем системы
				p1 = particles[w];
				dt = t_global - p1.t;
				p1.x[0] = p1.x[0] + p1.vx * dt;
				p1.y[0] = p1.y[0] + p1.vy * dt;
				p1.z[0] = p1.z[0] + p1.vz * dt;

				// Cуммируем значения координат центров узлов за первые m2 и за 
				// последние m2 соударений
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

				// Запоминаем все координаты всех частиц
				x[m][i] = p1.x[0];
				y[m][i] = p1.y[0];
				z[m][i] = p1.z[0];
			}
		}

		// Находим средние значения координат центров узлов за первые m2 и за 
		// последние m2 соударений
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

		// Начинаем поиск линии, по которой будем определять поворот 
		// осей координат для исследуемого слоя

		// Для начала находим самую левую нижнюю частицу в слое, из числа тех,
		// которые были отобраны для рассмотрения (не ближе чем 6 радиусов от
		// периодических граничных условий):
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

					// Ищем ближайшую к данной частице частицу на расстоянии не дальше 3.2 радиусов
					// так, чтобы линия, пересекающая центры их узлов,
					// образовывала острый угол с осью OY.
					if ((sqrt(dy*dy + dz*dz) < 3.2*particle_R) && (dy > 0.0 && dz > 0.0))
					{
						alfa = atan2(dy, dz);

						// Если такая частица найдена и она лежит на той же прямой, что
						// и все уже отобранные для этой линии ранее частицы, то записываем ее в список:
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

					// Ищем частицу на расстоянии не дальше 3.2 радиусов
					// так, чтобы линия, пересекающая центры узлов,
					// образовывала острый угол с осью OY.
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

		// Берем среднюю частицу в списке частиц вдоль одной из осей симметрии слоя
		// и находим для этой частицы другую ось симметрии, составляя список всех частиц,
		// центры которых находятся вдоль этой линии:
		min_m = list_of_particles_alfa[particles_count_alfa / 2];
		repeat = 1;
		while (repeat > 0)
		{
			repeat = -1;
			for (int m = 0; repeat < 0 && m < particles_count; m++) {
				if (m != min_m) {
					dy = middle_y1[min_m] - middle_y1[m];
					dz = middle_z1[min_m] - middle_z1[m];

					// Ищем частицу на расстоянии не дальше 3.2 радиусов
					// так, чтобы линия, пересекающая центры узлов,
					// образовывала острый угол с осью OY.
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

		// Получаем угол поворота линии методом наименьших квадратов
		// для всех частиц в массиве list_of_particles_alfa (это одна из
		// осей кристалла в данном слое, эта ось образует острый угол с осью OY
		// и проходит через 10+ узлов в данном слое):
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

		// K_alfa1 - это коэффициент в уравнении прямой (z = ky + b), рассчитанный на основании
		// положения частиц в данном слое на протяжении первых m2 соударений:
		K_alfa1 = (particles_count_alfa*sum_yz - sum_y*sum_z) / (particles_count_alfa*sum_zz - sum_z*sum_z);
		// K_alfa1 - это коэффициент в уравнении прямой (z = ky + b), рассчитанный на основании
		// положения частиц в данном слое на протяжении последних m2 соударений:
		K_alfa2 = (particles_count_alfa*sum_yz2 - sum_y2*sum_z2) / (particles_count_alfa*sum_zz2 - sum_z2*sum_z2);

		alfa = atan(K_alfa1);
		alfa2 = atan(K_alfa2);

		B_alfa1 = (sum_y - K_alfa1*sum_z) / particles_count_alfa;
		B_alfa2 = (sum_y - K_alfa2*sum_z) / particles_count_alfa;

		sum_y = 0.0; sum_z = 0.0; sum_yz = 0.0; sum_zz = 0.0;
		sum_y2 = 0.0; sum_z2 = 0.0; sum_yz2 = 0.0; sum_zz2 = 0.0;

		// Получаем угол поворота второй линии методом наименьших квадратов
		// (данная линия не совпадает с первой линией, но она пересекается с ней
		// в точке, совпадающей с центром одного из узлов, через которые проходит
		// первая линия):
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

		double delta_alfa = (alfa2 - alfa) / (m1 - m2); // Находим скорость вращения слоя (в градусах за время в среднем между
														// двумя измерениями)
		fprintf(debug_file, "\n /////////////////////////////////////////////////////\n");
		fprintf(debug_file, "\n delta alfa = %.5le\n", delta_alfa);

		// Находим пересечение двух найденных осей кристалла в данном слое
		// чтобы это сделать - решаем уравнение пересечения двух прямых, образованных
		// выбранными на предыдущих шагах частиц:
		// y1, z1 - координаты точки пересечения этих линий за первые 200 соударений,
		// y2, z2 - координаты точки пересечения линий за последние 200 соударений.
		double y1 = K_alfa1*(B_betta1 - B_alfa1) / (K_alfa1 - K_betta1) + B_alfa1;
		double y2 = K_alfa2*(B_betta2 - B_alfa2) / (K_alfa2 - K_betta2) + B_alfa2;
		double z1 = (B_betta1 - B_alfa1) / (K_alfa1 - K_betta1);
		double z2 = (B_betta2 - B_alfa2) / (K_alfa2 - K_betta2);

		fprintf(debug_file, "\n delta Y = %.5le \n delta Z = %.5le", y2 - y1, z2 - z1);

		fprintf(debug_file, "\n Y1 = %.5le, Z1 = %.5le \n", y1, z1);
		min_m = list_of_particles_alfa[particles_count_alfa / 2];
		fprintf(debug_file, "\n min_m Y = %.5le, Z = %.5le \n", middle_y1[min_m], middle_z1[min_m]);

		// Рассчитываем смещение всего слоя кристалла:
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

		// Проходим по всем запомненным нами координатам с m2/2 соударения
		// на частицу до m1-m2/2 соударения, то есть берем для рассмотрения m1-m2 соударений:
		for (int i = m2 / 2; i < m1 - m2 / 2; i++) {
			// Для каждой частицы в данном слое (если частица прошла все условия отбора):
			for (int j = 0; j < particles_count; j++) {

				// Запрашиваем номер частиц в списке частиц данного слоя:
				w = list_of_particles[j];

				// Рассчитываем номер текущей итерации чтобы сделать поправку на поворот и
				// смещение слоя кристалла:
				int f = i - m2 / 2;

				// Рассчитываем координаты частицы с учетом смещения слоя и среднего положения
				// центра данной частицы (то есть рассчитываем текущее смещение частицы на данном шаге от
				// среднего положения центра узла, который образован этой частицей):
				double x_mol = x[j][i] - (middle_x1[j] + middle_x2[j]) / 2.0;
				double y_mol0 = (middle_y1[j] - y[j][i] + delta_y*f)*cos(alfa) + (middle_z1[j] - z[j][i] + delta_z*f)*sin(alfa);
				double z_mol0 = -(middle_y1[j] - y[j][i] + delta_y*f)*sin(alfa) + (middle_z1[j] - z[j][i] + delta_z*f)*cos(alfa);

				double y_mol = y_mol0*cos(delta_alfa*f) + z_mol0*sin(delta_alfa*f);
				double z_mol = -y_mol0*sin(delta_alfa*f) + z_mol0*cos(delta_alfa*f);

				// Записываем номер частицы и рассчитанное смещение в файл:
				fprintf(f1_file, "%d %.5le %.5le %.5le\n", w, x_mol, y_mol, z_mol);
			}
		}
	}

	// Закрываем все файлы, которые мы редактировали:
	fclose(f1_file);
	fclose(debug_file);

	// Удаляем из памяти все динамические массивы:
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
Функция сжатия системы, позволяет сжимать или расширять систему до заданной плотности
с заданным максимальным шагом по плотности.

Аргументы:
compress_to_etta - итоговая плотность
delta_max - максимальный разрешённый шаг по L
KK1 - количество соударений на одну частицу в системе между
маленькими сжатиями по плотности
KK2 - количество сдвигов стенки, после которых необходимо провести
дополнительно KK3 соударений
KK3 - количество соударение на одну частицу в системе каждые KK2
сдвигов стенки
type - тип сжатия:
0 - пододвигать только левую стенку
1 - пододвигать только правую стенку
2 - пододвигать одновременно обе стенки на одинаковое расстояние
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
		// Находим идеальный радиус для частиц, решая уравнение
		// методом Кардано:
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

	// задаём плотность с точностью в 12 знаков
	while (abs(etta - compress_to_etta) > 1.0e-12) {
		if (type < 3) {
			if (etta < compress_to_etta) {
				max = -100.0;
				min = 1.1e+10;

				// Смещаем все частицы в текущее время системы чтобы синхронизовать
				// все частицы и затем находим координаты частиц, наиболее близко
				// расположенных от идеальных стенок, чтобы знать на сколько можно
				// пододвинуть стенку не коснувшись частиц.
				// Так же определяем минимальное растояние на котором частицы находятся
				// друг от друга, если изменение плотности системы производится
				// увеличением радиуса частиц.

				t_global = get_maximum_particle_time();

				for (long int i = 0; i < NP; ++i) {
					dt = t_global - particles[i].t;
					x = particles[i].x[0] + particles[i].vx * dt;

					if (DIFFUSE_RIGHT_WALL == true && x > L / 2 && particles[i].vx > 0)
						x = particles[i].x_max;

					if (x < min) min = x;
					if (x > max) max = x;
				}


				// Рассчитываем на сколько близко частицы находятся к первой
				// и второй идеальным стенкам

				min = min - particle_R;
				max = L - max - particle_R;

				// выбираем наименьшее и этих расстояний и сохраняем в deltaL
				deltaL = min;
				if (deltaL > max) deltaL = max;

				// сжимаем не впритык к частицам и не слишком быстро
				deltaL = deltaL / 1.1;
				if (deltaL < 0.1e-12) deltaL = 0.01e-30; // не двигаем стенку если частицы близко
				if (deltaL > delta_max) deltaL = delta_max;  // сдвигаем стенку не больше чем на delta_L

				dx = deltaL / 2.0;  // рассчитываем смещение для всех частиц


									// Если следущее смещение стенки сожмёт систему больше
									// чем требуется (т.е. относительная плотность частиц в системе etta
									// будет больше / меньше той, которая была указана а аргументах),
									// то необходимо задать точное значение L, чтобы плотность совпала
									// с требуемой плотностью.

				if (L - deltaL < L_ideal) {
					deltaL = 0.0;
					L = L_ideal;
					dx = 0.0;
				}
			}
			else {

				//Если систему необходимо расширить, то просто берём
				//наибольший разрешённый шаг для изменения коорднаты стенки
				//и сдвигаем идеальную стенку или две стенки.

				deltaL = L - L_ideal;
				if (deltaL < -delta_max) deltaL = -delta_max; // шаг расширения системы
				dx = deltaL / 2.0;
			}

			if (type != 1) {

				// Двигаем все частицы влево, таким образом пододвигая только левую стенку.
				// расстояние частиц до правой стенки не изменится, т.к. после перемещения частиц
				// влево мы сдвигаем обе стенки на то же расстояние - в итоге мы сдвинем все
				// частицы влево на deltaL/2.0.

				for (long int w = 0; w < NP; ++w) {
					for (int j = 0; j < 4; j++) {
						particles[w].x[j] -= dx;
					}
				}

				if (type == 0) {
					deltaL = dx;
				}
			}

			// изменяем систему - происходит мгновенное изменение координат двух стенок
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

				// Следим чтобы не пересечь стенки
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

		//Сохраняем состояние системы и снова загружаем его пересчитав новые
		//параметры и проведя необходимую инициализацию
		save("tmp");
		load_seed("tmp");

		// После каждого смещения стенки делаем указанное количество соударений
		// на частицу в системе
		for (long int i = 0; i < KK1; i++) {
			step();
		}

		compression_steps_done += 1;
		if (compression_steps_done % KK2 == 0) {

			for (long int i = 0; i < KK3; i++) {
				step();

				// очистить строку
				printf("\r%76c\r", ' ');

				printf("Relaxation: %ld / %ld, current etta=%.15le ", i, KK3, etta);
			}
		}

		// рассчитываем плотность после сжатия
		etta = (4.0 * PI * particle_R3 * NP) / (3.0 * A * A * (L - particle_Rx2));

		// очистить строку
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
		// очистить строку
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
Функция проведения эксперимента по указанному в
текстовом файле описанию.

Аргументы:
file_name - имя файла с описанием шагов эксперимента для программы.
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
		srand(time(NULL));  // перемешиваем генератор случайных чисел

							// Если необходимо создать новый посев
		if (str_command.compare("new") == 0) {
			double etta;   // средняя плотность системы
			long int particles_count;

			command_file >> particles_count;
			command_file >> etta;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

			new_seed(particles_count, etta);

			print_system_parameters();
		}
		// Если необходимо загрузить состояние системы из файла
		if (str_command.compare("load") == 0) {
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

			load_seed(parameter);

			printf("\n System was successfully loaded from file '%s' \n", parameter);

			print_system_parameters();
		}
		// Если необходимо рассчитать динамику системы в течении
		// указанного количества соударений
		if (str_command.compare("step") == 0) {
			command_file >> steps;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

			printf(" INFO: Step start. \n");

			start = clock();

			long int f[100];
			for (i = 0; i < 100; i++) f[i] = 0;
			double t_max;
			double dt, x;

			for (i = 0; i < steps; ++i) {
				step();

				// очистить строку
				printf("\r%76c\r", ' ');

				printf("%ld / %ld", long(i + 1), steps);

				// если мы пишем файл с историей событий,
				// то очищать его каждые 1000 соударений на частицу,
				// чтобы он не превышал допустимых размеров.
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
			Рассчитываем момент импульса системы Lx и выводим его на экран.
			*/
			double Lx = 0.0;
			for (i = 0; i < NP; i++) {
				Lx += particles[i].y[0] * particles[i].vz - particles[i].z[0] * particles[i].vy;
			}
			printf("\n Lx = %.15le \n", Lx);
			check_particles();
		}
		// если необходимо собрать данные по профилю плотности в системе
		if (str_command.compare("image") == 0) {
			command_file >> steps; // число соударений на одну частицу за время измерений
			command_file >> i; // точность. число точек графика на один радиус частицы
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

			image(steps, i, parameter);
		}
		// если необходимо собрать данные по "разрезу" в системе
		if (str_command.compare("profile") == 0) {
			int dots_for_each_particle;
			double x1, x2;
			command_file >> x1;
			command_file >> x2;
			command_file >> dots_for_each_particle;
			command_file >> steps;
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

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
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

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
		// если необходимо сохранить состояние системы
		if (str_command.compare("save") == 0) {
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

			save(parameter);
			printf("\n INFO: coordinates of particles saved to '%s' \n", parameter);
		}
		if (str_command.compare("function_g") == 0) {
			double x1, x2;
			command_file >> steps;
			command_file >> x1;
			command_file >> x2;
			command_file >> parameter;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

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
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

													   //count_of_nearest_particles(steps, x1, x2, parameter);
			printf("\n Function count_of_nearest_particles finished \n");
		}
		if (str_command.compare("measure_pressure") == 0) {
			command_file >> steps;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

			measure_pressure(steps);
		}
		if (str_command.compare("enable_autosave_for_compress") == 0) {
			temp_save = true;
			command_file >> temp_save_file;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки
		}
		if (str_command.compare("disable_autosave_for_compress") == 0) {
			temp_save = false;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки
		}
		// если необходимо изменить плотность системы
		if (str_command.compare("compress") == 0) {
			command_file.getline(parameter, 255, '\n');  // завершить считывание строки
			printf("\n Sorry, 'compress' command was deprecated, we need to use \n");
			printf("'compress_two_walls'(short form: 'compresst') instead. You can also use ");
			printf("'compress_left_wall'(short form: 'compressl') \n or 'compress_right_wall'");
			printf("(short form: 'compressr') commands to change system density.\n");
			exit(1);
		}
		// если необходимо изменить плотность системы
		// сжимать, сдвигая две идеальные стенки
		if ((str_command.compare("compress_two_walls") == 0) ||
			(str_command.compare("compresst") == 0)) {
			long int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // требуемая плотность
			command_file >> delta_L;  // минимальное допустимое значение перемещения стенки
			command_file >> KK1;  // количество соударений после каждого шага сжатия
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки
			compress(etta, delta_L, KK1, KK2, KK3, 2);

			print_system_parameters();
		}
		// если необходимо изменить плотность системы
		// двигать только "левую" идеальную стеку, x = 0
		if ((str_command.compare("compress_left_wall") == 0) ||
			(str_command.compare("compressl") == 0)) {
			long int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // требуемая плотность
			command_file >> delta_L;  // минимальное допустимое значение изменения плотности
			command_file >> KK1;  // количество соударений после каждого шага сжатия
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки
			compress(etta, delta_L, KK1, KK2, KK3, 0);

			print_system_parameters();
		}
		// если необходимо изменить плотность системы
		// двигать только "правую" идеальную стеку, x = L
		if ((str_command.compare("compress_right_wall") == 0) ||
			(str_command.compare("compressr") == 0)) {
			long int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // требуемая плотность
			command_file >> delta_L;  // минимальное допустимое значение изменения плотности
			command_file >> KK1;  // количество соударений после каждого шага сжатия
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки
			compress(etta, delta_L, KK1, KK2, KK3, 1);

			print_system_parameters();
		}
		// если необходимо изменить плотность системы
		// сжимать, изменяя радиус частиц
		if (str_command.compare("compress_R") == 0) {
			long int KK1, KK2, KK3;
			double etta, delta_L;
			command_file >> etta; // требуемая плотность
			command_file >> delta_L;  // минимальное допустимое значение изменения радиуса
			command_file >> KK1;  // количество соударений после каждого шага сжатия
			command_file >> KK2;
			command_file >> KK3;
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки
			compress(etta, delta_L, KK1, KK2, KK3, 3);

			print_system_parameters();
		}
		if (str_command.compare("diffuse_right_wall_on") == 0 ||
			str_command.compare("diff_on") == 0) {
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

			DIFFUSE_RIGHT_WALL = true;

			save("tmp.txt");
			load_seed("tmp.txt");
		}
		if (str_command.compare("diffuse_right_wall_off") == 0 ||
			str_command.compare("diff_off") == 0) {
			command_file.getline(tmp_str, 255, '\n');  // завершить считывание строки

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
Стартовая функция для программы, считать файл "program.txt"
и начать выполнять эксперимент, описанный в этом файле.
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