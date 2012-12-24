/*

  не трогать функции работы с очередью событий, не оптимизировать.
    они отлажены и должны быть именно такими.

TO DO:
 1. virtual particles

    ИМПОРТИРУЕМЫЕ МОДУЛИ
*/
#include <time.h>
#include <fstream>
#include <math.h>
/*
    ОБЪЯВЛЕНИЕ ПЕРЕМЕНЫХ
*/

// число ячееек, на которые разбивается объём 
short K, K2;

// число частиц во всём объёме
#define NP 6976

// число пи - для расчёта средней плотности системы
#define PI 3.14159265358979

// параметры объёма, задаются в load()
double A, A2, dA, L, dL;

// указатель на последний элемент в очереди событий.
short last;

// объект "событие"
typedef struct Event_ 
{
	double t;
	short im, jm;
} Event;

// очередь событий - оптимально 8192 элемента
Event time_queue[8192];		// или больше?

// объект "частица"
typedef struct particle_ {
	double x, y, z, vx, vy, vz, t, dt;
	short x_box, y_box, z_box, ti, box_i;
} particle;

// массив частиц
particle particles[2*NP];

// клетка. Объём системы разделём на множествво клеток,
// каждая клетка содержит в себе несколько виртуальных частиц
typedef struct Box_
{
	double x1, x2, y1, y2, z1, z2;
	short particles[11];
	short end;
} Box;

// массив клеток для всего объёма
Box boxes_yz[16][16][47];

/*
 функция подъема элемента по очереди событий к началу очереди
 параметры:
 i - позиция, на которой находится элемент в данный момент
 t - время до наступления данного события
*/
int get_up(short i, double &t)
{
	short int j = i >> 1;
	while (i > 1)
	{
		if (i % 2 != 0 && t < time_queue[i-1].t)
		{
			particles[time_queue[i-1].im].ti = i;
			if (time_queue[i-1].jm >= 0) particles[time_queue[i-1].jm].ti = i;
			time_queue[i] = time_queue[i-1];
			i--;
		}
		if (time_queue[j].t > t)
		{
			particles[time_queue[j].im].ti = i;
			if (time_queue[j].jm >= 0) particles[time_queue[j].jm].ti = i;
			time_queue[i] = time_queue[j];
			i = j; j >>= 1;
		} else return i;
	}
	return 1;
}

/*
 функция добавления события в очередь событий
 параметры:
 e - новое событие (см. подробнее структуру типа Event)
*/
void add_event(short &i, short &j)
{
	double t = particles[i].t + particles[i].dt;

	particles[i].ti = get_up(last, t);

	if (j >= 0) 
	{
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
void delete_event(short i)
{
	short j = i << 1;
	while (j < last)
	{
		if (i % 2 == 0 && time_queue[i+1].t < time_queue[j].t)
		{
			particles[time_queue[i+1].im].ti = i;
			if (time_queue[i+1].jm >= 0) particles[time_queue[i+1].jm].ti = i;
			time_queue[i] = time_queue[i+1];
			i++; j = i << 1;
		}
		else
		{
			particles[time_queue[j].im].ti = i;
			if (time_queue[j].jm >= 0) particles[time_queue[j].jm].ti = i;
			time_queue[i] = time_queue[j];
			i = j; j = i << 1;
		}
	}

	if (i < last - 1 && i % 2 == 0)
	{
		particles[time_queue[i+1].im].ti = i;
		if (time_queue[i+1].jm >= 0) particles[time_queue[i+1].jm].ti = i;
		time_queue[i] = time_queue[i+1];
		i++;
	}

	if (i < last-1)
	{
		j = get_up(i, time_queue[last-1].t);
		particles[time_queue[last-1].im].ti = j;
		if (time_queue[last-1].jm >= 0) particles[time_queue[last-1].jm].ti = j;
		time_queue[j] = time_queue[last-1];
	}
	
	last--;
}


void create_virt_particle(short &i, double &Y, double &Z)
{
    double x, y, z, t1, t2;
    int k = 0;
    int w = -1;
    particle p = particles[i];

    t1 = - Y / p.vy;
    t2 = - Z / p.vz;
    if (t1 > t2) t1 = t2;
    x = p.x - p.vx*t1;
    while (true)
    {
        if (x < 1.0) x = -x;
        else
        if (x > L-1.0) x = x-L+1.0;
        else
        break;
        ++k;
    }
    if (k % 2 == 1) { p.vx = -p.vx; p.vy = -p.vy; p.vz = -p.vz; }

    y = (p.vy/p.vx)*(x-p.x) + p.y;
    z = (p.vz/p.vx)*(x-p.x) + p.z;

    particles[i].x_box = short(x / dL + 1);
    particles[i].y_box = short(y / dA);
    particles[i].z_box = short(z / dA);

    for (s = 0; s <= boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].end; ++s)
	{
		n = boxes_yz[q2][w2][r].particles[s];
		temp = p.t - particles[n].t;

		dx = particles[n].x + particles[n].vx*temp - x;
		dy = particles[n].y + particles[n].vy*temp - y;
		dz = particles[n].z + particles[n].vz*temp - z;

        if (dx*dx + dy*dy + dz*dz < 4.0) w = n;
	}

    if (w == -1) //	вставляем частицу
    {
        particles[i].box_i = ++boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].end;
		boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].particles[particles[i].box_i] = i;
    }

                   // иначе: ищем частицу для столкновения и сталкиваем их.
}


/*
 функция расчёта близжаёшего события для частицы
 параметры:
 i - номер частицы
*/
void retime(short &i) 
{
	particle p1 = particles[i];
	short jm, bx, by, bz;
	double dt, dt2;
	Box box = boxes_yz[p1.y_box][p1.z_box][p1.x_box];

	if (p1.vx < 0.0) 
	{
		dt2 = (box.x1 - p1.x) / p1.vx;
		jm = -2;
		if (p1.x_box == 1) 
		{
			dt2 = (box.x1 + 1.0 - p1.x) / p1.vx;
			jm = -1;
		}
	}
	else 
	{
		dt2 = (box.x2 - p1.x) / p1.vx;
		jm = -4;
		if (p1.x_box == K2-1) 
		{
			dt2 = (box.x2 - 1.0 - p1.x) / p1.vx;
			jm = -1;
		}
	}

	if (p1.vy < 0.0) 
	{
		dt = (box.y1 - p1.y) / p1.vy;
		if (dt < dt2) { dt2 = dt; jm = -5; }
	}
	else 
	{
		dt = (box.y2 - p1.y) / p1.vy;
		if (dt < dt2) { dt2 = dt; jm = -6; }
	}

	if (p1.vz < 0.0) 
	{
		dt = (box.z1 - p1.z) / p1.vz;
		if (dt < dt2) { dt2 = dt; jm = -7; }
	}
	else 
	{
		dt = (box.z2 - p1.z) / p1.vz;
		if (dt < dt2) { dt2 = dt; jm = -8; }
	}

	double temp, dx, dy, dz, dvx, dvy, dvz, d, dv, bij;
	short s, n, r, q, w, q2, w2, i_dy, i_dz;
	for (r = p1.x_box - 1; r < p1.x_box + 2; ++r)
		for (q = p1.y_box - 1; q < p1.y_box + 2; ++q)
			for (w = p1.z_box - 1; w < p1.z_box + 2; ++w) 
			{
				q2 = q; w2 = w;
				if (q==-1) q2 = K-1;
				if (q==K) q2 = 0;
				if (w==-1) w2 = K-1;
				if (w==K) w2 = 0;

				for (s = 0; s <= boxes_yz[q2][w2][r].end; ++s)
				{
					n = boxes_yz[q2][w2][r].particles[s];
					temp = p1.t - particles[n].t;
					dvx = particles[n].vx - p1.vx;
					dvy = particles[n].vy - p1.vy;
					dvz = particles[n].vz - p1.vz;
					dx = particles[n].x + particles[n].vx*temp - p1.x;
					dy = particles[n].y + particles[n].vy*temp - p1.y;
					dz = particles[n].z + particles[n].vz*temp - p1.z;

					bij = dx*dvx + dy*dvy + dz*dvz;
					if (bij < 0.0)
					{
						dv = dvx*dvx + dvy*dvy + dvz*dvz;
						d = bij*bij + dv*(4.0 - dx*dx - dy*dy - dz*dz);
						if (d > 0.0)
						{
							dt = -(sqrt(d)+bij)/dv;
							temp += dt;

							if ((dt < dt2) && (temp < particles[n].dt) && (dt > 0.0))  
							{
								dt2 = dt; jm = n;
							}
						}
					}
				}
			}

	if (jm >= 0) 
	{
		dt = p1.t - particles[jm].t;
		particles[jm].t = p1.t; 
		particles[jm].dt = dt2;
		particles[jm].x += particles[jm].vx * dt; 
		particles[jm].y += particles[jm].vy * dt; 
		particles[jm].z += particles[jm].vz * dt;

		if (time_queue[particles[jm].ti].jm >= 0)
		{
			if (time_queue[particles[jm].ti].im == jm) 
				time_queue[particles[jm].ti].im = time_queue[particles[jm].ti].jm;
			time_queue[particles[jm].ti].jm = -100;
		} else delete_event(particles[jm].ti);
	}

	particles[i].dt = dt2;
	add_event(i, jm);
}

/*
 функция загрузки данных о системе из файла.
 загружаются параметры объёма, координаты и скорости всех частиц
*/
void load_seed() {
	double a1, a2, a3, a4, x, y, z;
	int i;
	FILE *loading_file = fopen("save_file.txt", "r");
	fscanf(loading_file, "%i\n", &i);
	fscanf(loading_file, "%le\n", &a1);
	A = a1 / 2.0; A2 = a1;
	fscanf(loading_file, "%le\n", &a1);
	L = a1 / 2.0;
	printf("N=%i  A=%f  L=%f \n", i, 2.0*A, 2.0*L);
	dA = A2 / K;
	dL = L*2.0 / (K2-1);
	y = z = -A;
	x = -L - dL;
	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
		{
			for (int w = 0; w <= K2; w++)
			{
				boxes_yz[i][j][w].x1 = x;
				boxes_yz[i][j][w].x2 = x+dL;
				boxes_yz[i][j][w].y1 = y;
				boxes_yz[i][j][w].y2 = y+dA;
				boxes_yz[i][j][w].z1 = z;
				boxes_yz[i][j][w].z2 = z+dA;
				boxes_yz[i][j][w].end = -1;
				x += dL;
			}
			x = -L - dL;
			z += dA;
		}
		y += dA;
		z = -A;
	}

	short end, x_box, y_box, z_box;
	for (short i = 0; i < NP; i++)
    {
		fscanf(loading_file, "%le %le %le %le\n", &a1, &a2, &a3, &a4);		
		particles[i].x = a1 - L; 
		particles[i].y = a2 - A;
		particles[i].z = a3 - A; 
		x_box = short(a1 / dL + 1);
		y_box = short(a2 / dA);
		z_box = short(a3 / dA);

		if (boxes_yz[y_box][z_box][x_box].x1 < particles[i].x &&
			boxes_yz[y_box][z_box][x_box].x2 > particles[i].x)
			particles[i].x_box = x_box;
		else
		{
			printf("A");
		}

		if (boxes_yz[y_box][z_box][x_box].y1 < particles[i].y &&
			boxes_yz[y_box][z_box][x_box].y2 > particles[i].y)
			particles[i].y_box = y_box;
		else
		{
			printf("B");
		}

		if (boxes_yz[y_box][z_box][x_box].z1 < particles[i].z &&
			boxes_yz[y_box][z_box][x_box].z2 > particles[i].z)
			particles[i].z_box = z_box;
		else
		{
			printf("C");
		}

        fscanf(loading_file, "%le %le %le\n", &a1, &a2, &a3);
		particles[i].vx = a1;
		particles[i].vy = a2;
		particles[i].vz = a3;
		particles[i].t = 0.0; particles[i].dt = 0.0; particles[i].ti = last;
		end = particles[i].box_i = ++boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].end;
		boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].particles[particles[i].box_i] = i;

		if (end > 7)
			printf("alarm");
		for (short t = 0; t < particles[i].box_i; ++t)
			for (short d = t+1; d <= particles[i].box_i; ++d)
				if (boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].particles[t] == boxes_yz[particles[i].y_box][particles[i].z_box][particles[i].x_box].particles[d])
					printf("alarm");
	}
	fclose(loading_file);
	double etta = (PI * NP) / (6.0 * A * A * (L - 1.0) );
	printf("etta=%.15le\n", etta);

	last = 1;
	time_queue[0].t = 0.0;
	for (short i = 1; i < 6000; ++i) time_queue[i].t = 1.0E+20;

	for (short i = 0; i < NP-1; ++i)
		for (short j = i+1; j < NP; ++j)
		{
			double x = particles[i].x - particles[j].x;
			double y = particles[i].y - particles[j].y;
			double z = particles[i].z - particles[j].z;
			double d = x*x + y*y + z*z;
			if (d < 3.99999)
			{
				printf("\n %.15le\n", d);
				printf("R");
			}
		}

	for (short i = 0; i < NP; ++i) retime(i);
}

/*
 функция изменения состояния частиц в соответвии с наступившим
 в системе событием.
 параметры:
 e - наступившее событие
*/
void reform(short &im, short &jm) 
{
	particle p1 = particles[im];
	double dx, dy, dz, q1, q2, z;
    short bx, by, bz;

	bx = p1.x_box; by = p1.y_box; bz = p1.z_box;

	if (jm >= 0)
	{
		particle p2 = particles[jm];
		p2.x += p2.vx * p2.dt; p2.y += p2.vy * p2.dt; p2.z += p2.vz * p2.dt;
		p2.t = p1.t; p2.dt = 0.0;
		dx = p1.x - p2.x;
		dy = p1.y - p2.y;
		dz = p1.z - p2.z;
		q1 = (dx*p1.vx + dy*p1.vy + dz*p1.vz) / 4.0;
		q2 = (dx*p2.vx + dy*p2.vy + dz*p2.vz) / 4.0;
		z = q2-q1; p1.vx += dx*z; p1.vy += dy*z; p1.vz += dz*z;
		z = q1-q2; p2.vx += dx*z; p2.vy += dy*z; p2.vz += dz*z;
		particles[jm] = p2;
	}
	else
	if (jm == -1)
	{
		p1.vx = -p1.vx;
	}
	else 
	{
		if (jm != -100)
		{
			short end = boxes_yz[by][bz][bx].end;
			boxes_yz[by][bz][bx].particles[p1.box_i] = boxes_yz[by][bz][bx].particles[end];
			particles[boxes_yz[by][bz][bx].particles[end]].box_i = p1.box_i;
			--boxes_yz[by][bz][bx].end;

			if (jm == -2)
			{
				--p1.x_box;
				p1.x = boxes_yz[by][bz][p1.x_box].x2;
			}
			if (jm == -4)
			{
				++p1.x_box;
				p1.x = boxes_yz[by][bz][p1.x_box].x1;
			}
			if (jm == -5)
			{
				--p1.y_box;
				if (by < 0) p1.y_box = K-1;
				p1.y = boxes_yz[p1.y_box][bz][bx].y2;
			}
			if (jm == -6)
			{
				++p1.y_box;
				if (p1.y_box == K) p1.y_box = 0;
				p1.y = boxes_yz[p1.y_box][bz][bx].y1;
			}
			if (jm == -7)
			{			
				--p1.z_box;
				if (p1.z_box < 0) p1.z_box = K-1;
				p1.z = boxes_yz[by][p1.z_box][bx].z2;
			}
			if (jm == -8)
			{
				++p1.z_box;
				if (p1.z_box == K) p1.z_box = 0;
				p1.z = boxes_yz[by][p1.z_box][bx].z1;
			}
			p1.box_i = ++boxes_yz[p1.y_box][p1.z_box][p1.x_box].end;
			boxes_yz[p1.y_box][p1.z_box][p1.x_box].particles[p1.box_i] = im;
		}
	}
	particles[im] = p1;
}

/*
 функция "шаг", основной цикл программы
*/
void step() 
{
	particle p1;
	int i;
	short im, jm, r=0;
	double time = 0.0;

	while (r < NP)
	{
		im = time_queue[1].im;
		jm = time_queue[1].jm;

		delete_event(1);

		p1 = particles[im];

		if (jm == -100) p1.dt = time - p1.t;
		p1.x += p1.vx * p1.dt;
		p1.y += p1.vy * p1.dt;
		p1.z += p1.vz * p1.dt;
		p1.t += p1.dt; p1.dt = 0.0; time = p1.t;
		particles[im] = p1;
		reform(im, jm);
		if (jm >= 0) { retime(jm); ++r; }
		retime(im);
	}

	particle *p = particles;
	for (i = 0; i < NP; ++i, ++p)
		(*p).t -= time;
	Event *t = time_queue; ++t;
	for (i = 1; i < last; ++i, ++t)
		(*t).t -= time;
	time = 0.0;
}

/*
 функция сохранения состояния системы в текстовый файл.
 перед сохранением производится синхронизация частиц по времени.
 параметры:
 f - имя файла, в который будет записана информация
*/
void save(char f[255])
{
	double x, y, z;

	FILE *save_file = fopen(f, "w+");
    fprintf(save_file, "%d\n", NP);
    fprintf(save_file, "%.15le\n", A*2.0);
    fprintf(save_file, "%.15le\n", L*2.0);
	for (short int i = 0; i < NP; ++i) {
		x = L + particles[i].x - particles[i].vx * particles[i].t;
		y = A + particles[i].y - particles[i].vy * particles[i].t;
		z = A + particles[i].z - particles[i].vz * particles[i].t;
		fprintf(save_file, "%.15le %.15le %.15le %.15le %.15le \n", x, y, z, particles[i].t, particles[i].dt);
		fprintf(save_file, "%.15le %.15le %.15le \n", particles[i].vx, particles[i].vy, particles[i].vz);
    }
	fclose(save_file);
}


/*
void init()
{
	using namespace std;
	clock_t start, end, result;
	char command[255], file_name[255];
	int i, steps;
	ifstream command_file("program.txt");

	K = 10;
	K2 = 44;

	while (!command_file.eof())
	{
		command_file.getline(command, 255, ' ');

		if (command[0] == 'K')
		{
			if (command[1] == 'X') { command_file >> K2; command_file.getline(command, 5, '\n'); }
			if (command[1] == 'Y') { command_file >> K; command_file.getline(command, 5, '\n'); }
		}
		
		if (command[0] == 'l' && command[1] == 'o' &&
			command[2] == 'a' && command[3] == 'd') { load_seed(); command_file.getline(command, 5, '\n'); }
 		if (command[0] == 's' && command[1] == 't') 
		{
			command_file >> steps;
			command_file.getline(command, 5, '\n');

			System::Console::Write("\nStep start " + steps + "\n");

			start = clock();
			
			for (i=0; i < steps; ++i) step();
			
			end = clock();
			result = end - start;

			System::Console::Write(" \7\a\n ======\n " + steps + " collisions per particle\n");
			System::Console::Write(" Total Time = " + (result/CLOCKS_PER_SEC) + " seconds.\n ======\n");
		}
		if (command[0] == 's' && command[1] == 'a')
		{
			command_file.getline(file_name, 255);
			save(file_name);
		}
	}
}
*/


void image()
{
	int k; int W = (int(L)+1)*40;
	printf("L .%f W %d\n", L, W);
	int img[2000];
	FILE *save_file = fopen("img.txt", "w+");

	printf("image started \n");
	for (short g = 0; g < 2000; ++g) img[g] = 0;
	for (short h = 0; h < 1000; ++h)
	{
		step();

		for (short i = 0; i < NP; ++i) 
        {
			k = int((L + particles[i].x - particles[i].vx * particles[i].t)*200.0)/10;
			++img[k];
		}
	}
	for (short f = 0; f < W; ++f)
	{
		fprintf(save_file, "%d\n", img[f]);
	}
	fclose(save_file);

	printf("image complete\n");
}


void f1()
{
	clock_t start, end, result;

	printf("10000 collissions per particle\n");

	start = clock();
	for (short i = 0; i < 10000; ++i) 
	{
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d/10000", i);
		step();
	}
	end = clock(); result = end - start;
	
	printf("\nComplete. Total Time = %d seconds\n", (int)(result/CLOCKS_PER_SEC) );

	save("save_file.tut");
    	load_seed();
	image();
	printf("image complete");
}

/*
 начальная точка входа.
 засекаем время, запускаем основной цикл и ждем окончания его выполнения
*/
int main()
{
	K = 12; K2 = 46;
	load_seed();

    printf("\n\n %d \n\n", K2);
 
    f1();
	
	return 0;
}
