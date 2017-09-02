//============================================================================
// Name        : fun.cpp
// Author      : l
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cstdio>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <ctime>

#include "lib3ds/file.h"
#include "lib3ds/io.h"
#include "lib3ds/mesh.h"

#include "BVH.h"
#include "Triangle.h"
#include "ComplexVector3.h"

using std::vector;

using namespace std;

//#define Nx 0.9
//#define Ny 0.9
//const float dx = 0.02 / 1;
//const float dy = 0.02 / 1;
//const double PI = 3.14159265358979323846;
const double PI = 3.14159265358979323846;
const double c = 299792458;
const CXF j = CXF(0, 1);
const int maxhitnumber = 3;

class RayIndex {
public:
	int ray_one;
	int ray_two;
	int ray_three;
	int ray_four;
	int ray_center;

};


//
//static const int xdnum = Nx / dx;
//static const int ydnum = Ny / dy;
//int rnum = xdnum * ydnum;
//int tnum = (xdnum - 1) * (ydnum - 1);

//int GeoOpt(const Vector3 &raypoint, const BVH &bvh, Vector3 &hitpoint,
//		Vector3 &k, Vector3 &E, float &t);
bool intersectPlane(const Vector3 &n, const Vector3 &p0, const Vector3 &l0,
		const Vector3 &l, double &t);
void ComplexField(const Vector3 &R, const double &kr, const double frequency, CXV3 &C);
double Vector_Angle_d(const Vector3 &v1, const Vector3 &v2);
void rotate_THETA_then_PHI(Vector3* R1, int n, double THETA, double PHI);
//Vector3 generate_Ant(const float &THETA, const float &PHI,
//		const float &distance, Vector3 *antr, Vector3 *antt);
void genrateAperture(const double&THETA, const double &PHI, const double &distance,
		const double cellSize, const int rnum,const int tnum,
		const double &Lx, const double &Ly, const int xdnum, const int ydnum, Vector3 *antr, Vector3 *antt,
		Vector3  &ant_N);
void equation_twenty_seven(const CXV3 &E, const double &k0, const double cellSize, const Vector3 retpoint, const Vector3 &theta_h,
		const Vector3 &phi_h, CXF &Atheta, CXF &Aphi);
double SafeAcos(double x);

void calEquationTwentySeven(const CXV3 *E, const double &k0,
		const double &cellSize, const int &tnum, const Vector3 &theta_h,
		const Vector3 &phi_h, const double &PHI, const double *IR_Angle, Vector3 *retpoint_antr,
		CXF &Atheta, CXF &Aphi);

void AllrayGeoOpt(const int &tnum, const Vector3 *tube, const BVH &bvh, const Vector3 &kin,
		const Vector3 &Ein, Vector3 *refpoint_tube, Vector3 *k_tube, Vector3 *E_tube,
		double *kdr_tube, int *hitnum_tube, int &OneAngleNumberOfHit);

int GeoOpt(const Vector3 &raypoint, const BVH &bvh, const Vector3 &kin, const Vector3 &Ein,
		Vector3 &hitpoint, Vector3 &kdir, Vector3 &Eout, double &t);

void getMatlabData(const string &POL,  const int &PointOfAngle, double *RCSTheta, double *RCSPhi);


void initializeEF(const string &POL, const Vector3 &kin, Vector3 &theta_h, Vector3 &phi_h,
		Vector3 &Ei);

//void AllEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in, Vector3 *k_tube,
//		const Vector3 *tube, const Vector3 *refpoint_tube, const double *kdr_tube,
//		const double f, const Vector3 *E_tube,
//		double *kdr_tuber, Vector3 *retpoint_tube, double *total_path_tube,
//		CXV3 *Ec_tube, double *IR_Angle);
void AllEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in, const Vector3 *k_tube, Vector3 *k_antr,
		RayIndex *index,const int xdnum,
		const Vector3 *tube, const Vector3 *refpoint_tube, const double *kdr_tube,
		const double f, const Vector3 *E_tube,
		double *kdr_tuber, Vector3 *retpoint_tube, double *total_path_tube,
		CXV3 *Ec_tube, double *IR_Angle);
void Raytube_Numbering(const int &tnum, const int &xdnum, const Vector3 *antr, const Vector3 *antt, RayIndex *index) ;

void getAngleData(const int &PointOfAngle, const int *OneAngleNumberOfHit);

//


int main(int argc, char *argv[]) {
	double f = 10e9;
	double lambda = c / f;

	double Lx = 0.2;
	double Ly = 0.2;
	double rayPerWavelenth = 64;
	double cellSize = lambda/rayPerWavelenth;
	double Pmin = -45;
	double Pmax = 45;
	double Pstep = 0.1;
	int PointOfAngle = (Pmax-Pmin)/Pstep + 1;

	int xdnum = Lx / cellSize;
	int ydnum = Ly / cellSize;
	double k0 = 2 * PI * f * pow(c, -1);



	int rnum = xdnum * ydnum;
	int tnum = (xdnum - 1) * (ydnum - 1);

	Vector3* antr = new Vector3[rnum]; //ray number
	Vector3* tube = new Vector3[tnum]; //raytube number
	Vector3* k_antr = new Vector3[rnum];
	Vector3* k_tube = new Vector3[tnum];
	Vector3* E_antr = new Vector3[rnum];
	Vector3* E_tube = new Vector3[tnum];
	Vector3* kt_in = new Vector3[tnum];

	int *hitnum_antr = new int[rnum];
	int *hitnum_tube = new int[tnum];
	double* kdr_antr = new double[rnum]; // Point of emission to the point of reflection distance
	double* kdr_tube = new double[tnum];  //Point of emission to the point of reflection distance
	double* kdr_antrr = new double[rnum]; //reflection point
	double* kdr_tuber = new double[tnum]; //reflection point to reflection plane distance
	double* total_path_tube = new double[tnum];
	double* IR_Angle = new double[tnum];


	CXV3* Ec_antr = new CXV3[rnum];
	CXV3* Ec_tube = new CXV3[tnum];
	double* RCSTheta = new double[PointOfAngle];
	double* RCSPhi =new double[PointOfAngle];
	Vector3* retpoint_antr = new Vector3[rnum]; //The Return point on incident plane
	Vector3* refpoint_antr = new Vector3[rnum]; //The source point on the object
	Vector3* retpoint_tube = new Vector3[tnum]; //The Return point on incident plane
	Vector3* refpoint_tube = new Vector3[tnum]; //The source point on the object
	CXF* Atheta = new CXF[tnum];
	CXF* Aphi = new CXF[tnum];
	CXF* Atheta_total = new CXF[1];
	CXF* Aphi_total = new CXF[1];
	int* OneAngleNumberOfHit = new int[PointOfAngle];


	bool *return_hit = new bool[tnum];

	memset(refpoint_antr, 0, rnum);
	memset(hitnum_tube, 0, tnum);
	memset(kt_in, 0, tnum);

	RayIndex *index;
	index = new RayIndex[tnum];

//	string input = "/home/user/cuda-workspace/lib3/dihedral m.3ds";
//	string input = "/home/user/cuda-workspace/lib3/dihedral 22pt5d.3ds";
//	string input = "/home/user/cuda-workspace/lib3/dihedral 45d.3ds";
//	string input = "/home/user/cuda-workspace/lib3/Untitled800.3ds";
	string input = "/home/user/cuda-workspace/lib3/plate1515.3ds";



	Lib3dsFile* fin = lib3ds_file_load(input.c_str());

	if (!fin) {
		cerr << "*****ERROR****\n Loading file failed" << input << endl;
		exit(EXIT_FAILURE);
	} else {
		cout << "Read ok" << endl;
	}
	Lib3dsMesh* mesh = fin->meshes;
	size_t num_face = mesh->faces;
//	size_t num_vert = mesh->points;

	cout<<"#Face ="<<num_face<<endl;
//	cout<<"#Vertex="<<num_vert<<endl;

	Lib3dsPoint* vert = mesh->pointL;
	Lib3dsFace* face = mesh->faceL;

//	for(size_t i=0; i<num_vert; i++){
//		cout<<vert[i].pos[2]<<endl;
//	}
//	exit (0);

	vector<Object*> objects;

	for (size_t i = 0; i < num_face; i++) {
		objects.push_back(
				new Triangle(
						Vector3(vert[face[i].points[0]].pos[0],
								vert[face[i].points[0]].pos[1],
								vert[face[i].points[0]].pos[2]),
						Vector3(vert[face[i].points[1]].pos[0],
								vert[face[i].points[1]].pos[1],
								vert[face[i].points[1]].pos[2]),
						Vector3(vert[face[i].points[2]].pos[0],
								vert[face[i].points[2]].pos[1],
								vert[face[i].points[2]].pos[2])));
	}
	BVH bvh(&objects);

//	float THETA = 90.0 * PI / 180;
//	float PHI = 44 * PI / 180;
	double distance = 0.600;
	double PHI;
	double THETA = -90 * PI / 180;
	string POL = "H";
	Vector3 phi_h;
	Vector3 theta_h;
	Vector3 Ei;


	for (int i =0  ; i < PointOfAngle; i++) {
		PHI = (-i*Pstep + Pmax) * PI / 180;
		genrateAperture(THETA, PHI, distance, cellSize, rnum, tnum, Lx,
				Ly, xdnum, ydnum, antr, tube, kt_in[0]);




		initializeEF(POL, kt_in[0], theta_h, phi_h, Ei);



		AllrayGeoOpt(tnum, tube ,bvh, kt_in[0], Ei, refpoint_tube,  k_tube,  E_tube, kdr_tube,
				hitnum_tube, OneAngleNumberOfHit[i]);


//		AllrayGeoOpt(rnum, antr ,bvh, kt_in[0], Ei, refpoint_antr,  k_antr,  E_antr, kdr_antr,
//				hitnum_antr);




		 AllEF(tnum, hitnum_tube, kt_in, k_tube, k_antr, index, xdnum,
				tube, refpoint_tube, kdr_tube,
				f, E_tube,
				kdr_tuber, retpoint_tube, total_path_tube,
				Ec_tube, IR_Angle);

//			for (int n = 0; n<tnum ; n++){
//				intersectPlane((-1 * kt_in[0]), tube[0], refpoint_tube[n],  k_tube[n], kdr_tuber[n]);
//				retpoint_tube[n] = refpoint_tube[n] + kdr_tuber[n] * k_tube[n];
//
//				intersectPlane((-1 * kt_in[0]), tube[0], refpoint_antr[n],  k_antr[n], kdr_antrr[n]);
//				retpoint_antr[n] = refpoint_antr[n] + kdr_antrr[n] * k_antr[n];
//
//				}


	 calEquationTwentySeven(Ec_tube, k0,
				cellSize, tnum, theta_h,
				phi_h, PHI, IR_Angle, retpoint_antr, Atheta_total[0], Aphi_total[0]);


		RCSTheta[i] = 10* log10(PI * 4 * abs(Atheta_total[0]) * abs(Atheta_total[0]));
		RCSPhi[i]= 10 * log10(PI * 4 * abs(Aphi_total[0]) * abs(Aphi_total[0]));

		cout << "RCS(phi) ="
				<< RCSTheta[i]
				<< " dBsm" << endl;
		cout << "RCS(theta) ="<<
				RCSPhi[i]<< " dBsm"
				<< endl;


	}
	getAngleData(PointOfAngle, OneAngleNumberOfHit);

	getMatlabData(POL, PointOfAngle, RCSTheta, RCSPhi);





	delete[] antr;
	delete[] tube;
	delete[] k_antr;
	delete[] k_tube;
	delete[] E_antr;
	delete[] E_tube;
	delete[] hitnum_antr;
	delete[] hitnum_tube;
	delete[] kdr_antr; //distance Point of emission to the point of reflection
	delete[] kdr_tube; //Point of emission to the point of reflection
	delete[] kdr_antrr; //reflection point
	delete[] kdr_tuber; //reflection point

	delete[] Ec_antr;
	delete[] Ec_tube;

	delete[] retpoint_antr; //The Return point on incident plane
	delete[] refpoint_antr; //The source point on the object
	delete[] retpoint_tube; //The Return point on incident plane
	delete[] refpoint_tube; //The source point on the object
	delete[] Atheta;
	delete[] Aphi;
	delete[] total_path_tube;
	delete[] return_hit;
	delete[] index;
	delete[] kt_in;
	delete[] Atheta_total;
	delete[] Aphi_total;
	delete[] RCSTheta;
	delete[] RCSPhi;
	delete[] IR_Angle;
	delete[] OneAngleNumberOfHit;

	cout << (double)clock() / CLOCKS_PER_SEC << " S";
	return 0;

}


int GeoOpt(const Vector3 &raypoint, const BVH &bvh, const Vector3 &kin, const Vector3 &Ein,
		Vector3 &hitpoint, Vector3 &kdir, Vector3 &Eout, double &t) {

	Vector3 phi_h, theta_h, Ei, Xc, Yc, Zc, theta_ci_h, phi_ci_h, theta_cr_h,
			phi_cr_h, Nz, Ed, m_h, Er;
	double theta_ci;
//	float phi_ci;
	int hitnumber = 0;
	Ei = Ein;
	t = 0;

	Ray ray(raypoint, kin);
	IntersectionInfo I;


	bool hit = bvh.getIntersection(ray, &I, false);

	while (hit) {
		hitnumber = hitnumber + 1;

//		cout << "hitpoint = " << I.hit.x << "," << I.hit.y << "," << I.hit.z<< endl;
//                   cout<<I.t<<"I.t"<<endl;
		Vector3 nor = I.object->getNormal(I);
//                 cout<<nor.x<<","<<nor.y<<","<<nor.z<<"nor"<<endl;

		//R1 formula (11)
		theta_ci = acos((-1 * ray.d) * nor);
//		phi_ci = 0;
//		cout<<"theta_ci = "<<theta_ci<<endl;

		if (theta_ci > 1e-15) {
			//R1 formula (6)
			m_h = (-1 * ray.d) ^ nor / sin(theta_ci);
			//		               cout<<m_h.x<<","<<m_h.y<<","<<m_h.z<<"m"<<endl;

			//R1 formula (10)
			Xc = m_h ^ nor;
			Zc = -1 * nor;
			Yc = -1 * m_h;

			//R1 formula (12)
			theta_ci_h = Xc * cos(theta_ci) - Zc * sin(theta_ci);
			phi_ci_h = Yc;
			//R1 formula (13)
			theta_cr_h = (-1 * Xc * cos(theta_ci)) - Zc * sin(theta_ci);
			phi_cr_h = Yc;

			// get a reflection vector
			//using https://math.stackexchange.com/questions/13261/how-to-get-a-reflection-vector
			kdir = ray.d - 2 * (ray.d * nor) * nor;
			//    cout<<r.x<<","<<r.y<<","<<r.z<<"=k "<<endl;

			ray = Ray(I.hit, kdir);

			//R1 formula (8) and (9)
			if (hitnumber == 1) {

				Ed = (Ei * phi_ci_h) * phi_ci_h
						+ (Ei * theta_ci_h) * theta_ci_h;
				Er = (-1) * (Ei * phi_ci_h) * phi_cr_h
						+ (-1) * (Ei * theta_ci_h) * theta_cr_h;
			} else {
				Ed = Er;
				Er = (-1) * (Ed * phi_ci_h) * phi_cr_h
						+ (-1) * (Ed * theta_ci_h) * theta_cr_h;

			}

		} else {
			kdir = (-1) * ray.d;
			if (hitnumber % 2 == 1) {
				Er = (-1) * Ei;
			} else {
				Er = Ei;
			}
			ray = Ray(I.hit, kdir);
		}

		Eout = Er;
		hitpoint = I.hit;
		t += I.t;
//		cout<<Ei.x<<","<<Ei.y<<","<<Ei.z<<" = ei"<<endl;
//		cout << "hitpoint = " << hitpoint.x << "," << hitpoint.y << ","
//				<< hitpoint.z << endl;

//		    cout<<Ed.x<<","<<Ed.y<<","<<Ed.z<<"ed"<<endl;
//		cout << Er.x << "," << Er.y << "," << Er.z << " =er " << endl;
//		cout << hitnumber << "hitnumber" << endl;

		bool hit = bvh.getIntersection(ray, &I, false);

		if ((!hit) || (hitnumber == maxhitnumber)) {
			break;
		}

	}
	return hitnumber;
}
//Ray-Plane Intersection
//Ref: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
bool intersectPlane(const Vector3 &n, const Vector3 &p0, const Vector3 &l0,
		const Vector3 &l, double &t) {
	// assuming vectors are all normalized
	double denom = n * l;
	if (denom > 1e-15) {
		Vector3 p0l0 = p0 - l0;
		t = p0l0 * n / denom;
		return (t >= 0);
	}

	return false;
}
void ComplexField(const Vector3 &R, const double &kr, const double frequency, CXV3 &C) {
	double multiply2 = 2.0;
	CXF expvalue = exp(j * (kr) / c * frequency * multiply2 * PI);
	C = CXV3(CXF(R.x, 0) * expvalue, CXF(R.y, 0) * expvalue,
			CXF(R.z, 0) * expvalue);
//	C = CXV3(CXF(R.x, 0) , CXF(R.y, 0) , CXF(R.z, 0) );

}

double Vector_Angle_d(const Vector3 &v1, const Vector3 &v2) {
	double l1 = length(v1);
	double l2 = length(v2);
	double angle = 0;
//	cout << v2.x << "," << v2.y << "," << v2.z << endl;

//	cout<<acos(-1) <<endl;

	angle = SafeAcos(v1 * v2 / (l1 * l2));
	angle = angle * 180 / PI;
//	cout << angle << endl;

	return angle;
}

double SafeAcos(double x) {
	if (x < -1.0)
		x = -1.0;
	else if (x > 1.0)
		x = 1.0;
	return acos(x);
}

void rotate_THETA_then_PHI(Vector3 *R1, int n, double THETA, double PHI) {
//	//Roate X Axis (THETA)
//	Vector3* R2 = new Vector3[n];
//	for (int i = 0; i < n; i++) {
//		R2[i].x = R1[i].x;
//		R2[i].y = cos(THETA) * R1[i].y - sin(THETA) * R1[i].z;
//		R2[i].z = sin(THETA) * R1[i].y + cos(THETA) * R1[i].z;
//	}
	//Roate Y Axis (THETA)
	Vector3* R2 = new Vector3[n];

	for (int i = 0; i < n; i++) {
		R2[i].x = cos(THETA) * R1[i].x - sin(THETA) * R1[i].z;
		R2[i].y = R1[i].y;
		R2[i].z = sin(THETA) * R1[i].x + cos(THETA) * R1[i].z;
	}

	//Roate Z Axis (PHI)
	for (int i = 0; i < n; i++) {
		R1[i].x = cos(PHI) * R2[i].x - sin(PHI) * R2[i].y;
		R1[i].y = sin(PHI) * R2[i].x + cos(PHI) * R2[i].y;
		R1[i].z = R2[i].z;
	}
	delete[] R2;
}

//Vector3 generate_Ant(const float &THETA, const float &PHI,
//		const float &distance, Vector3 *antr, Vector3 *antt) {
//	Vector3 antc, mov, cros_1, cros_2, ant_N, ant_Nd;
//	Vector3 origin = Vector3(0, 0, 0);
//	float f = 15e9;
//	float lambda = c / f;
////	float Nx = 0.5;
////	float Ny = 0.5;
//	float dx = lambda / 64;
//	float dy = lambda / 64;
//	int xdnum = Nx / dx;
//	int ydnum = Ny / dy;
//	int rnum = xdnum * ydnum;
//	int tnum = (xdnum - 1) * (ydnum - 1);
//	for (int i = 0; i < rnum; i++) {
//		antr[i].x = 0 + (i % xdnum) * dx;
//		antr[i].y = 0 + i / xdnum * dy;
//		antr[i].z = 0;
//	}
////
////	for (int i = 0; i < rnum; i++) {
////		  cout<<antr[i].x<<","<<antr[i].y<<","<<antr[i].z<<"antr"<<endl;
////	}
//
//	for (int i = 0; i < tnum; i++) {
//		antt[i].x = (0 + (i % (xdnum - 1)) * dx + dx / 2);
//		antt[i].y = (0 + i / (xdnum - 1) * dy + dy / 2);
//		antt[i].z = 0;
//	}
//
//	antc = (antr[0] + antr[rnum - 1]) / 2.0;
//	mov = antc - origin;
//
//	for (int i = 0; i < rnum; i++) {
//		antr[i] = antr[i] - mov;
//	}
//
//	for (int i = 0; i < tnum; i++) {
//		antt[i] = antt[i] - mov;
//	}
//	antc = antc - mov;
////	cout<<antc.x<<","<<antc.y<<","<<antc.z<<endl;
////	cout<<mov.x<<","<<mov.y<<","<<mov.z<<endl;
//	rotate_THETA_then_PHI(antr, rnum, THETA, PHI);
//	rotate_THETA_then_PHI(antt, tnum, THETA, PHI);
//
//	cros_1 = antr[1] - antr[0];
//	cros_2 = antt[tnum - 1] - antr[0];
//
//	ant_N = cros_2 ^ cros_1;
//	ant_N = normalize(ant_N);
//	ant_Nd = distance * ant_N;
//
//	for (int i = 0; i < rnum; i++) {
//		antr[i] = antr[i] - ant_Nd;
//	}
//
//	for (int i = 0; i < tnum; i++) {
//
//		antt[i] = antt[i] - ant_Nd;
//	}
//
//	antc = antc - ant_Nd;
//	return ant_N;
//
//}

void AllrayGeoOpt(const int &tnum, const Vector3 *tube, const BVH &bvh, const Vector3 &kin,
		const Vector3 &Ein, Vector3 *refpoint_tube, Vector3 *k_tube, Vector3 *E_tube,
		double *kdr_tube, int *hitnum_tube, int &OneAngleNumberOfHit){
	for (int n = 0; n < tnum; n++) {
//			k_tube[n] =
//			E_tube[n] = Ei;
		//		cout<<n<<endl;
		hitnum_tube[n] = GeoOpt(tube[n], bvh,  kin, Ein, refpoint_tube[n], k_tube[n],
				E_tube[n], kdr_tube[n]);
		OneAngleNumberOfHit+=hitnum_tube[n];

	}

}
void equation_twenty_seven(const CXV3 &E, const double &k0, const double cellSize,
		const Vector3 retpoint, const Vector3 &theta_h,
		const Vector3 &phi_h, CXF &Atheta, CXF &Aphi) {



//	float Sx = 0;
//	float Sy = 0;
//	float Xi = retpoint * phi_h;
//	float Yi = retpoint * theta_h;
//	cout<<Xi<<","<<Yi<<endl;
//	cout << "E=" << E.x << "," << E.y << "," << E.z << endl;
//	CXF expvalue = exp(j * k0 * (Sx * Xi * +Sy * Yi));
	CXF expvalue = CXF(1,0);
	CXF Ex;
	Ex = E * phi_h;
//	if ( Ex!=CXF(0,0)){
//	cout<<"Ex ="<<Ex<<endl;
//	}
	CXF Ey;
	Ey = E * theta_h;
//	if ( Ey!=CXF(0,0)){
// 	cout<<"Ey ="<<Ey<<endl;
//	}
//	cout<<"Ey ="<<Ey<<endl;

//	float Phi_p; //Phi for paper
	double Ii = 1.0;

	if (Ex == CXF(0, 0) && Ey == CXF(0, 0)) {

		Atheta = CXF(0, 0);
		Aphi = CXF(0, 0);
	} else {

		Atheta = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii);
//			cout<<"Atheta= "<<Atheta<<endl;
		Aphi = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii);
//			cout<<"Aphi= "<<Aphi<<endl;

	}



}

void calEquationTwentySeven(const CXV3 *E, const double &k0,
		const double &cellSize, const int &tnum, const Vector3 &theta_h,
		const Vector3 &phi_h, const double &PHI, const double *IR_Angle, Vector3 *retpoint_antr,
		CXF &Atheta, CXF &Aphi){
	CXF* SingleRayAtheta = new CXF[tnum];
	CXF* SingleRayAphi = new CXF[tnum];
	CXF expvalue = CXF(1,0);
	CXF Ex, Ey;
	double Ii = 1.0;

	Atheta = CXF(0, 0);
	Aphi = CXF(0, 0);
	for(int i = 0 ; i < tnum ; i++){

			Ex = E[i] * phi_h;
		//	if ( Ex!=CXF(0,0)){
		//	cout<<"Ex ="<<Ex<<endl;
		//	}
			Ey = E[i] * theta_h;
		//	if ( Ey!=CXF(0,0)){
		// 	cout<<"Ey ="<<Ey<<endl;
		//	}
		//	cout<<"Ey ="<<Ey<<endl;

		//	float Phi_p; //Phi for paper

			if (Ex == CXF(0, 0) && Ey == CXF(0, 0)) {

				SingleRayAtheta[i] = CXF(0, 0);
				SingleRayAphi[i] = CXF(0, 0);
			} else {

		//		SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii)/cos((180-IR_Angle[i])*PI/180);
				SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii)/cos(PHI)/cos(PHI);
		//		SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii)/cos(PHI);

		//	   SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii)*(1+tan(PHI)*tan(2*PHI));
		//	   SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii)*(1+tan(PHI)*tan(2*PHI));
		//	   SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii);

		//			cout<<"Atheta= "<<Atheta<<endl;
		//		SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii)/cos((180-IR_Angle[i])*PI/180);
				SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii)/cos(PHI)/cos(PHI);
		//		SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii)/cos(PHI);

		//		SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii)/cos(IR_Angle[i]*PI/180);
		//		SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii)*(1+tan(PHI)*tan(2*PHI));
		//		SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii);



		//			cout<<"Aphi= "<<Aphi<<endl;

			}



			Atheta += SingleRayAtheta[i];
			Aphi += SingleRayAphi[i];
	}
	delete[] SingleRayAtheta;
	delete[] SingleRayAphi;


}
void getAngleData(const int &PointOfAngle, const int *OneAngleNumberOfHit){
	ofstream angledata("angledata.out");
	for (int i = 0; i < PointOfAngle; i++) {
		 angledata <<  OneAngleNumberOfHit[i] <<endl;

		ifstream ifile("angledata.out");
	}
	angledata.close();
}

void getMatlabData(const string &POL,  const int &PointOfAngle, double *RCSTheta, double *RCSPhi){

		if(POL =="V"){

		ofstream file("ThetapolRcsTheta.out");
		for (int i = 0; i < PointOfAngle; i++) {
			file << RCSTheta[i] <<endl;
	//		file << RCSPhi[i] << endl;
			ifstream ifile("ThetapolRcsTheta.out");
		}
		file.close();

		ofstream file1("ThetapolRcsPhi.out");
		for (int i = 0; i < PointOfAngle; i++) {
	//		file << RCSTheta[i] <<endl;
			file1 << RCSPhi[i] << endl;
			ifstream ifile1("ThetaPolRcsPhi.out");
		}
		file1.close();
		}else{
			ofstream file("PhipolRcsTheta.out");
			for (int i = 0; i < PointOfAngle; i++) {
				file << RCSTheta[i] <<endl;
		//		file << RCSPhi[i] << endl;

				ifstream ifile("PhipolRcsTheta.out");
			}
			file.close();

			ofstream file1("PhipolRcsPhi.out");
			for (int i = 0; i < PointOfAngle ; i++) {
		//		file << RCSTheta[i] <<endl;
				file1 << RCSPhi[i] << endl;
				ifstream ifile1("PhipolRcsPhi.out");
			}
			file1.close();
		}
}

void genrateAperture(const double&THETA, const double &PHI, const double &distance,
		const double cellSize, const int rnum,const int tnum,
		const double &Lx, const double &Ly, const int xdnum, const int ydnum, Vector3 *antr, Vector3 *antt,
		Vector3  &ant_N){
	Vector3 antc, cros_1, cros_2, ant_Nd;

	for (int i = 0; i < rnum; i++) {
		antr[i].x = 0 + (i % xdnum) * cellSize;
		antr[i].y = 0 + i / xdnum * cellSize;
		antr[i].z = 0;
	}

	for (int i = 0; i < tnum; i++) {
		antt[i].x = (0 + (i % (xdnum - 1)) * cellSize + cellSize / 2);
		antt[i].y = (0 + i / (xdnum - 1) * cellSize + cellSize / 2);
		antt[i].z = 0;
	}
//
//	for (int i = 0; i < rnum; i++) {
//		  cout<<antr[i].x<<","<<antr[i].y<<","<<antr[i].z<<"antr"<<endl;
//	}

	antc = (antr[0] + antr[rnum - 1]) / 2.0;


	for (int i = 0; i < rnum; i++) {
		antr[i] = antr[i] - antc;
	}


	for (int i = 0; i < tnum; i++) {
		antt[i] = antt[i] - antc;
	}

	rotate_THETA_then_PHI(antr, rnum, THETA, PHI);
	rotate_THETA_then_PHI(antt, tnum, THETA, PHI);


	cros_1 = antr[1] - antr[0];
	cros_2 = antr[rnum] - antr[0];

	ant_N = cros_2 ^ cros_1;
	ant_N = normalize(ant_N);
	ant_Nd = distance * ant_N;

	for (int i = 0; i < rnum; i++) {
		antr[i] = antr[i] - ant_Nd;
	}

	for (int i = 0; i < tnum; i++) {

		antt[i] = antt[i] - ant_Nd;
	}


//	antc = antc - ant_Nd;

}

void initializeEF(const string &POL, const Vector3 &kin, Vector3 &theta_h, Vector3 &phi_h,
		Vector3 &Ei){

	Vector3 Nz = Vector3(0, 0, 1);

	 phi_h = kin ^ Nz;
	 theta_h = kin ^ phi_h;
//		cout << "phi_h = " << phi_h.x << "," << phi_h.y << "," << phi_h.z
//				<< endl;
//		cout << "theta_h = " << theta_h.x << "," << theta_h.y << ","
//				<< theta_h.z << endl;

	if (POL == "V") {
		Ei = theta_h; //v ->Y'
	} else {
		Ei = phi_h;   //h ->X'
	}
	//	 Ei = (Ei * theta_h) * theta_h + (Ei * phi_h) * phi_h;

}

void AllEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in, const Vector3 *k_tube, Vector3 *k_antr,
		RayIndex *index,const int xdnum,
		const Vector3 *tube, const Vector3 *refpoint_tube, const double *kdr_tube,
		const double f, const Vector3 *E_tube,
		double *kdr_tuber, Vector3 *retpoint_tube, double *total_path_tube,
		CXV3 *Ec_tube, double *IR_Angle){

	int a, b;
//	for (int i = 0; i < tnum; i++) {
//		a = i / (xdnum - 1);
//		b = i % (xdnum - 1);
//
//		index->ray_center = i;
//		(*index).ray_one = xdnum * a + b;
//		(*index).ray_two = (*index).ray_one + 1;
//		(*index).ray_three = (*index).ray_two + xdnum;
//		(*index).ray_four = (*index).ray_one + xdnum;
//
//	}
	for (int n = 0; n < tnum; n++) {
		a = n / (xdnum - 1);
		b = n % (xdnum - 1);

		index->ray_center = n;
		(*index).ray_one = xdnum * a + b;
		(*index).ray_two = (*index).ray_one + 1;
		(*index).ray_three = (*index).ray_two + xdnum;
		(*index).ray_four = (*index).ray_one + xdnum;


		if(hitnum_tube[n]!=0){
			IR_Angle[n] = Vector_Angle_d(kt_in[0], k_tube[n]);
//				cout<<IR_Angle[n]<<endl;

		}else{
			IR_Angle[n] = -1;
		}

//
//		if (IR_Angle[n] > 45 && k_tube[index->ray_center].x==k_antr[index->ray_four].x &&
//				k_tube[index->ray_center].y==k_antr[index->ray_four].y &&
//				k_tube[index->ray_center].z==k_antr[index->ray_four].z &&
//				 k_tube[index->ray_center].x==k_antr[index->ray_one].x &&
//				 k_tube[index->ray_center].y==k_antr[index->ray_one].y &&
//				 k_tube[index->ray_center].z==k_antr[index->ray_one].z &&
//				 k_tube[index->ray_center].x==k_antr[index->ray_three].x &&
//				 k_tube[index->ray_center].y==k_antr[index->ray_three].y &&
//				 k_tube[index->ray_center].z==k_antr[index->ray_three].z &&
//				 k_tube[index->ray_center].x==k_antr[index->ray_two].x &&
//				 k_tube[index->ray_center].y==k_antr[index->ray_two].y &&
//				 k_tube[index->ray_center].z==k_antr[index->ray_two].z
//
//		)
		if (IR_Angle[n] > 45 ){
			intersectPlane((-1 * kt_in[0]), tube[0],
					refpoint_tube[n], (-1 * kt_in[0]), kdr_tuber[n]);
//			intersectPlane((-1 * kt_in[0]), tube[0],
//					refpoint_tube[n], k_tube[n], kdr_tuber[n]);


			retpoint_tube[n] = refpoint_tube[n] + kdr_tuber[n] * k_tube[n];
//			total_path_tube[n] = kdr_tuber[n] + kdr_tube[n];
			total_path_tube[n] = 2*kdr_tube[n];

//							cout << total_path_tube[n] << "totalpath" << endl;

			ComplexField(E_tube[n], total_path_tube[n], f, Ec_tube[n]);
		} else {
			Ec_tube[n] = CXV3((CXF(0, 0)), CXF(0, 0), CXF(0, 0));
		}
		//		cout << "rethit = " << return_hit[n] << endl;
	}

//	for (int n=0; n<tnum; n++){
//		if(hitnum_tube[n]!=0){
//			total_path_tube[n] = 2* kdr_tube[n];
//			ComplexField(E_tube[n], total_path_tube[n], f, Ec_tube[n]);
//		}else{
//			Ec_tube[n] = CXV3((CXF(0, 0)), CXF(0, 0), CXF(0, 0));
//		}
//
//	}
}

//void Raytube_Numbering(const int &tnum, const int &xdnum, const Vector3 *antr, const Vector3 *antt, RayIndex *index) {
//
////	RayIndex *index;
////	index = new RayIndex[tnum];
//
//	int m, n;
//	for (int i = 0; i < tnum; i++) {
//		m = i / (xdnum - 1);
//		n = i % (xdnum - 1);
//
//		index->ray_center = i;
//		(*index).ray_one = xdnum * m + n;
//		(*index).ray_two = (*index).ray_one + 1;
//		(*index).ray_three = (*index).ray_two + xdnum;
//		(*index).ray_four = (*index).ray_one + xdnum;
//
////
////    cout << index->ray_center<<","\
////         << (*index).ray_one << ","\
////         << (*index).ray_two << ","\
////         << (*index).ray_three << ","\
////         << (*index).ray_four << endl;
////
////
////
////  cout << antt[index->ray_center].x << ", " << antt[index->ray_center].y << "," << antt[index->ray_center].z  << "tube" << endl;
////  cout << antr[(*index).ray_one].x << ", " << antr[(*index).ray_one].y << "," << antr[(*index).ray_one].z  << "antr_one" << endl;
////  cout << antr[(*index).ray_two].x << ", " << antr[(*index).ray_two].y << "," << antr[(*index).ray_two].z  << "antr_two" << endl;
////  cout << antr[(*index).ray_three].x << ", " << antr[(*index).ray_three].y << "," << antr[(*index).ray_three].z  << "antr_three" << endl;
////  cout << antr[(*index).ray_four].x << ", " << antr[(*index).ray_four].y << "," << antr[(*index).ray_four].z  << "antr_four" << endl;
//
//	}
//}




