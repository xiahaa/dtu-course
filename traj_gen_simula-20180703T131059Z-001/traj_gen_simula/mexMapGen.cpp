#include "mex.h" 
#include "perlinnoise.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

typedef struct _PointXYZ
{
	float x;
	float y;
	float z;
}PointXYZ;

void randomMapGenerate(int sizeX, int sizeY, int sizeZ, double scale,
	int seed, int _ObsNum, double _w_l, double _w_h,
	std::vector<PointXYZ> &points);
	
void perlin3D(int sizeX, int sizeY, int sizeZ, double scale,
	int seed, double complexity, double fill, int fractal, double attenuation,
	std::vector<PointXYZ> &points);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{ 
	int seed = mxGetScalar(prhs[0]);
	int sizeX = mxGetScalar(prhs[1]);
	int sizeY = mxGetScalar(prhs[2]);
	int sizeZ = mxGetScalar(prhs[3]);
	double scale = mxGetScalar(prhs[4]);
	int type = mxGetScalar(prhs[5]);

	scale = 1 / scale;
	sizeX = sizeX * scale;
	sizeY = sizeY * scale;
	sizeZ = sizeZ * scale;
	std::vector<PointXYZ> points;

	if (type == 2)
	{
		double _w_l = mxGetScalar(prhs[6]);
		double _w_h = mxGetScalar(prhs[7]);
		int    _ObsNum = mxGetScalar(prhs[8]);

		mexPrintf("sx %d, sy %d, sz %d, resolution %f, type %d, w_l:%f, w_h:%f, obsnum:%d/n",
			sizeX,sizeY,sizeZ,scale,type,_w_l,_w_h,_ObsNum); 

		randomMapGenerate(sizeX, sizeY, sizeZ, scale, seed, _ObsNum, _w_l, _w_h, points);
	}
	else if (type == 1)
	{
		double complexity;
		double fill;
		int    fractal;
		double attenuation;
		if (nrhs == 10)
		{
			complexity = mxGetScalar(prhs[6]);
			fill = mxGetScalar(prhs[7]);
			fractal = mxGetScalar(prhs[8]);
			attenuation = mxGetScalar(prhs[9]);
		}
		else
		{
			complexity = 0.142857;
			fill = 0.38;
			fractal = 1;
			attenuation = 0.5;
		}
		mexPrintf("params %d complexity %f, fill %f, fractal %d, attenuation %f\n", nrhs, complexity,
			fill, fractal, attenuation);
		perlin3D(sizeX, sizeY, sizeZ, scale, seed, complexity, fill, fractal, attenuation, points);
	}
	
	mexPrintf("hello,world!/n"); 
	
	int N = points.size();
	plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL);
	plhs[1]=mxCreateDoubleMatrix(N,1,mxREAL);
	plhs[2]=mxCreateDoubleMatrix(N,1,mxREAL);
	
	double *ox, *oy, *oz;
	ox = mxGetPr(plhs[0]);
	oy = mxGetPr(plhs[1]);
	oz = mxGetPr(plhs[2]);
	for (int i = 0; i < N; i++)
	{
		ox[i] = points[i].x;
		oy[i] = points[i].y;
		oz[i] = points[i].z;		
	}
} 

void randomMapGenerate(int sizeX, int sizeY, int sizeZ, double scale,
	int seed, int _ObsNum, double _w_l, double _w_h,
	std::vector<PointXYZ> &points)
{
	std::default_random_engine eng(seed);

	double _resolution = 1 / scale;

	double _x_l = -sizeX / (2 * scale);
	double _x_h = sizeX / (2 * scale);
	double _y_l = -sizeY / (2 * scale);
	double _y_h = sizeY / (2 * scale);
	double _h_l = 0;
	double _h_h = sizeZ / scale;

	std::uniform_real_distribution<double> rand_x;
	std::uniform_real_distribution<double> rand_y;
	std::uniform_real_distribution<double> rand_w;
	std::uniform_real_distribution<double> rand_h;

	rand_x = std::uniform_real_distribution<double>(_x_l, _x_h);
	rand_y = std::uniform_real_distribution<double>(_y_l, _y_h);
	rand_w = std::uniform_real_distribution<double>(_w_l, _w_h);
	rand_h = std::uniform_real_distribution<double>(_h_l, _h_h);

	PointXYZ pt_random;
	points.clear();
	for (int i = 0; i < _ObsNum; i++)
	{
		double x, y;
		x = rand_x(eng);
		y = rand_y(eng);

		double w, h;
		w = rand_w(eng);
		h = rand_h(eng);

		int widNum = ceil(w / _resolution);
		int heiNum = ceil(h / _resolution);

		int rl, rh, sl, sh;
		rl = -widNum / 2;
		rh = widNum / 2;
		sl = -widNum / 2;
		sh = widNum / 2;

		for (int r = rl; r < rh; r++)
			for (int s = sl; s < sh; s++)
			{
				for (int t = 0; t < heiNum; t++)
				{
					if ((r - rl) * (r - rh + 1) * (s - sl) * (s - sh + 1) * t * (t - heiNum + 1) == 0)
					{
						pt_random.x = x + r * _resolution;
						pt_random.y = y + s * _resolution;
						pt_random.z = t * _resolution;
						points.push_back(pt_random);
					}
				}
			}
	}
}

void perlin3D(int sizeX, int sizeY, int sizeZ, double scale,
	int seed, double complexity, double fill, int fractal, double attenuation,
	std::vector<PointXYZ> &points)
{
	int width = sizeX * sizeY * sizeZ;
	int height = 1;
	points.resize(width * height);

	PerlinNoise noise(seed);

	std::vector<double>* v = new std::vector<double>;
	v->reserve(width);
	for (int i = 0; i < sizeX; ++i)
	{
		for (int j = 0; j < sizeY; ++j)
		{
			for (int k = 0; k < sizeZ; ++k)
			{
				double tnoise = 0;
				for (int it = 1; it <= fractal; ++it)
				{
					int    dfv = pow(2, it);
					double ta  = attenuation / it;
					tnoise += ta *
							noise.noise(dfv * i * complexity, dfv * j * complexity, dfv * k * complexity);
				}
				v->push_back(tnoise);
			}
		}
	}
	std::sort(v->begin(), v->end());
	int    tpos = width * (1 - fill);
	double tmp  = v->at(tpos);
	mexPrintf("threshold: %lf", tmp);

	int pos = 0;
	for (int i = 0; i < sizeX; ++i)
	{
		for (int j = 0; j < sizeY; ++j)
		{
			for (int k = 0; k < sizeZ; ++k)
			{
				double tnoise = 0;
				for (int it = 1; it <= fractal; ++it)
				{
					int    dfv = pow(2, it);
					double ta  = attenuation / it;
					tnoise += ta *
						noise.noise(dfv * i * complexity, dfv * j * complexity, dfv * k * complexity);
				}
				if (tnoise > tmp)
				{
					points[pos].x = i / scale - sizeX / (2 * scale);
					points[pos].y = j / scale - sizeY / (2 * scale);
					points[pos].z = k / scale;
					pos++;
				}
			}
		}
	}
	width = pos;
	points.resize(width * height);
}