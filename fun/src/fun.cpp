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
const float PI = 3.14159265358979323846;
const float c = 299792458;
const CXF j = CXF(0, 1);


const int maxhitnumber = 3;
//
//static const int xdnum = Nx / dx;
//static const int ydnum = Ny / dy;
//int rnum = xdnum * ydnum;
//int tnum = (xdnum - 1) * (ydnum - 1);

int GeoOpt(const Vector3 &raypoint, const BVH &bvh, Vector3 &hitpoint,
		Vector3 &k, Vector3 &E, float &t);
bool intersectPlane(const Vector3 &n, const Vector3 &p0, const Vector3 &l0,
		const Vector3 &l, float &t);
void ComplexField(const Vector3 &R, const float &kr, const float frequency, CXV3 &C);
float Vector_Angle_d(const Vector3 &v1, const Vector3 &v2);
void rotate_THETA_then_PHI(Vector3* R1, int n, float THETA, float PHI);
//Vector3 generate_Ant(const float &THETA, const float &PHI,
//		const float &distance, Vector3 *antr, Vector3 *antt);
void GenrateAperture(const float&THETA, const float &PHI, const float &distance,
		const float cellSize, const int rnum,const int tnum,
		const float &Lx, const float &Ly, const int xdnum, const int ydnum, Vector3 *antr, Vector3 *antt,
		Vector3  &ant_N);
//void Raytube_Numbering(Vector3 *antr, Vector3 *antt);
void equation_twenty_seven(const CXV3 &E, const float &k0, const float cellSize, const Vector3 retpoint, const Vector3 &theta_h,
		const Vector3 &phi_h, CXF &Atheta, CXF &Aphi);
double SafeAcos(double x);

class RayIndex {
public:
	int ray_one;
	int ray_two;
	int ray_three;
	int ray_four;
	int ray_center;

};








int main(int argc, char *argv[]) {
	float f = 15e9;
	float lambda = c / f;
//	float Nx = 0.5;
//	float Ny = 0.5;
	float Lx = 0.9;
	float Ly = 0.9;
	float rayPerWavelenth = 64;
	float cellSize = lambda/rayPerWavelenth;

//	float dx = lambda / 64;
//	float dy = lambda / 64;
	int xdnum = Lx / cellSize;
	int ydnum = Ly / cellSize;
	float k0 = 2 * PI * f * pow(c, -1);
	Vector3 Nz = Vector3(0, 0, 1);



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
	float* kdr_antr = new float[rnum]; // Point of emission to the point of reflection distance
	float* kdr_tube = new float[tnum];  //Point of emission to the point of reflection distance
	float* kdr_antrr = new float[rnum]; //reflection point
	float* kdr_tuber = new float[tnum]; //reflection point to reflection plane distance
	float* total_path_tube = new float[tnum];
	float* IR_Angle = new float[tnum];


//	CXV3* Ec_antr = new CXV3[rnum];
	CXV3* Ec_tube = new CXV3[tnum];
	double* RCSTheta = new double[tnum];
	double* RCSPhi =new double[tnum];
	Vector3* retpoint_antr = new Vector3[rnum]; //The Return point on incident plane
	Vector3* refpoint_antr = new Vector3[rnum]; //The source point on the object
	Vector3* retpoint_tube = new Vector3[tnum]; //The Return point on incident plane
	Vector3* refpoint_tube = new Vector3[tnum]; //The source point on the object
	CXF* Atheta = new CXF[tnum];
	CXF* Aphi = new CXF[tnum];
	CXF* Atheta_total = new CXF[1];
	CXF* Aphi_total = new CXF[1];


	bool *return_hit = new bool[tnum];

	memset(refpoint_antr, 0, rnum);
	memset(hitnum_tube, 0, tnum);
	memset(kt_in, 0, tnum);

	RayIndex *index;
	index = new RayIndex[tnum];

//	string input = "/home/user/cuda-workspace/lib3/dihedral m.3ds";
//	string input = "/home/user/cuda-workspace/lib3/dihedral 22pt5d.3ds";
//	string input = "/home/user/cuda-workspace/lib3/dihedral 45d.3ds";
	string input = "/home/user/cuda-workspace/lib3/Untitled800.3ds";


	Lib3dsFile* fin = lib3ds_file_load(input.c_str());

	if (!fin) {
		cerr << "*****ERROR****\n Loading file failed" << input << endl;
		exit(EXIT_FAILURE);
	} else {
		cout << "Read ok" << endl;
	}
	Lib3dsMesh* mesh = fin->meshes;
	size_t num_face = mesh->faces;
	size_t num_vert = mesh->points;

	cout<<"#Face ="<<num_face<<endl;
	cout<<"#Vertex="<<num_vert<<endl;

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
	float distance = 0.5500;
	float PHI;
	float THETA = -90 * PI / 180;
	string POL = "H";


	for (int i = 450 ; i <= 450; i++) {
		PHI = (-i*0.1 + 45) * PI / 180;

		//	Vector3 k;

//		kt_in[0] = generate_Ant(THETA, PHI, distance, antr, tube);

		GenrateAperture(THETA, PHI, distance, cellSize, rnum, tnum,Lx,
				Ly, xdnum, ydnum, antr, tube, kt_in[0]);

		//  for (int i = 0; i < rnum; i++) {
		//
		//    cout << antr[i].x << "," << antr[i].y << "," << antr[i].z << "antr"
		//        << endl;
		//  }
		//	for (int i = 0; i < tnum; i++) {
		//
		//		cout << tube[i].x << "," << tube[i].y << "," << tube[i].z << "tube"
		//				<< endl;
		//	}

		//	Vector3 k = Vector3(-1, 0, 0);
//		cout << "ki = " << " [ " << kt_in[0].x << " , " << kt_in[0].y << " , "
//				<< kt_in[0].z << " ] " << endl;
		Vector3 Ei;

		Vector3 phi_h = kt_in[0] ^ Nz;
		Vector3 theta_h = kt_in[0] ^ phi_h;
//		cout << "phi_h = " << phi_h.x << "," << phi_h.y << "," << phi_h.z
//				<< endl;
//		cout << "theta_h = " << theta_h.x << "," << theta_h.y << ","
//				<< theta_h.z << endl;


		if (POL == "V") {
			Ei = theta_h; //v ->Y'
		} else {
			Ei = phi_h;   //h ->X'
		}
//		cout << "Ei = " << " [ " << Ei.x << " , " << Ei.y << " , " << Ei.z
//				<< " ] " << endl;
//
//		Ei = (Ei * theta_h) * theta_h + (Ei * phi_h) * phi_h;
//		cout << "Ei = " << " [ " << Ei.x << " , " << Ei.y << " , " << Ei.z
//				<< " ] " << endl;



		//	for (int n = 0; n < rnum; n++) {
		//		k_antr[n] = k;
		//		E_antr[n] = Ei;
		//		hitnum_antr[n] = GeoOpt(antr[n], bvh, refpoint_antr[n], k_antr[n],
		//				E_antr[n], kdr_antr[n]);
		//	}

		for (int n = 0; n < tnum; n++) {
			k_tube[n] = kt_in[0];
			E_tube[n] = Ei;
			//		cout<<n<<endl;
			hitnum_tube[n] = GeoOpt(tube[n], bvh, refpoint_tube[n], k_tube[n],
					E_tube[n], kdr_tube[n]);
		}
//		for (int n = 0; n < tnum; n++) {
//			if (hitnum_tube[n]!=0){
//				cout << kdr_tube[n] << "=kdr" << endl;
//				cout << kdr_tuber[n] << "=kdrr" << endl;
//
//				cout << k_tube[n].x << "," << k_tube[n].y << "," << k_tube[n].z << "=k"
//										<< endl;
//				cout << E_tube[n].x << "," << E_tube[n].y << "," << E_tube[n].z << "=Er"<< endl;
//
//			}
//		}


		//	Raytube_Numbering(antr, tube);

		//	cout << "===============tube========================" << endl;


		//	cout << tnum << endl;
		for (int n = 0; n < tnum; n++) {
			if(hitnum_tube[n]!=0){
				IR_Angle[n] = Vector_Angle_d(kt_in[0], k_tube[n]);
//				cout<<IR_Angle<<endl;
			}else{
				IR_Angle[n] = -1;
			}

			if (IR_Angle[n] > 90) {
				return_hit[n] = intersectPlane((-1 * kt_in[0]), tube[0],
						refpoint_tube[n], k_tube[n], kdr_tuber[n]);

				retpoint_tube[n] = refpoint_tube[n] + kdr_tuber[n] * k_tube[n];
				total_path_tube[n] = kdr_tuber[n] + kdr_tube[n];
//							cout << total_path_tube[n] << "totalpath" << endl;

				ComplexField(E_tube[n], total_path_tube[n], f, Ec_tube[n]);
			} else {
				Ec_tube[n] = CXV3((CXF(0, 0)), CXF(0, 0), CXF(0, 0));
			}

			//		cout << "rethit = " << return_hit[n] << endl;

		}

//		for (int n = 0; n < tnum; n++) {
//			if (hitnum_tube[n]!=0){
//	            cout<<"PHI = "<<PHI*180/PI<<endl;
//	            cout<<"incident point ="<<"["<<tube[n].x<<","<<tube[n].y<<","<<tube[n].z<<"]"<<endl;
//	            cout<<"IR_Angle = "<<IR_Angle[n]<<endl;
//				cout<<"Number of hit = "<<hitnum_tube[n]<<endl;
//                cout <<  "total path = "<< total_path_tube[n] << endl;
//                cout <<"return point form target = "<<"["<<refpoint_tube[n].x<<","<<refpoint_tube[n].y
//                		<<","<<refpoint_tube[n].z<<"]"<<endl;
//                cout<<"return point = "<<"["<< retpoint_tube[n].x<<","<<retpoint_tube[n].y<<
//                		","<<retpoint_tube[n].z<<"]"<<endl;
//
//
//				cout << "sorce to target path = "<<kdr_tube[n] << endl;
//				cout <<"target return path = "<< kdr_tuber[n]  << endl;
//				cout <<"return k = "<<" [ "<< k_tube[n].x << " , " << k_tube[n].y << " , "
//						<< k_tube[n].z <<" ] "<< endl;
//				exit(0);
//
//			}
//
//		}


//		for (int n = 0; n < tnum; n++) {

			//		cout << kdr_tube[n] << "=kdr" << endl;
//					cout << kdr_tuber[n] << "=kdrr" << endl;
			//
			//		cout << k_tube[n].x << "," << k_tube[n].y << "," << k_tube[n].z << "=k"
			//				<< endl;
//					cout << Ec_tube[n].x << "," << Ec_tube[n].y << "," << Ec_tube[n].z << "=Ec"<< endl;

//		}

//		for (int n = 0; n < tnum; n++) {
//					cout << retpoint_tube[n].x << "," << retpoint_tube[n].y << ","
//							<< retpoint_tube[n].z << endl;

			//		cout << Ec_tube[n] << " = E" << endl;

//		}


		for (int n = 0; n < tnum; n++) {


			equation_twenty_seven(Ec_tube[n], k0, cellSize, retpoint_tube[n],
					theta_h, phi_h, Atheta[n], Aphi[n]);

			float denom2 = 2;
			Atheta[n] = j * k0 * Atheta[n] / PI / denom2;
			Aphi[n] = j * k0 * Aphi[n] / PI / denom2;
//					cout << "Aphi=";
//					cout << Aphi[n] << endl;
		}
		Atheta_total[0] = CXF(0, 0);
		Aphi_total[0] = CXF(0, 0);

		for (int n = 0; n < tnum; n++) {
			Atheta_total[0] += Atheta[n];
			Aphi_total[0] += Aphi[n];

		}

//		cout << "complex Atheta = " << Atheta_total[0] << endl;
//		cout << "complex Aphi = " << Aphi_total[0] << endl;

		//	cout<<PI*4*abs(Atheta_total[0])<<endl;
		cout << "RCS(phi) ="
				<< 10 * log10(PI * 4 * abs(Aphi_total[0]) * abs(Aphi_total[0]))
				<< " dBsm" << endl;
		cout << "RCS(theta) ="<<
				 10* log10(PI * 4 * abs(Atheta_total[0])* abs(Atheta_total[0])) << " dBsm"
				<< endl;
		RCSTheta[i] = 10* log10(PI * 4 * abs(Atheta_total[0]) * abs(Atheta_total[0]));
		RCSPhi[i]=10 * log10(PI * 4 * abs(Aphi_total[0]) * abs(Aphi_total[0]));

//		float sigma = 8 * PI * 0.3 * 0.3 * 0.3 * 0.3 / lambda / lambda;
//		sigma = 4 * PI * 0.8 * 0.8 * 0.8 * 0.8 / lambda / lambda;

//		sigma = 10 * log10(sigma);
//			cout << "sigma = " << sigma << " dBsm" << endl;
//			cout<<"lambda = "<<lambda<< "m"<<endl;

	}
//	cout<<POL<<endl;


	if(POL == "V"){

	ofstream file("ThetapolRcsTheta.out");
	for (int i = 0; i <= 900; i++) {
		file << RCSTheta[i] <<endl;
//		file << RCSPhi[i] << endl;

		ifstream ifile("ThetapolRcsTheta.out");
	}
	file.close();

	ofstream file1("ThetapolRcsPhi.out");
	for (int i = 0; i <= 900; i++) {
//		file << RCSTheta[i] <<endl;
		file1 << RCSPhi[i] << endl;
		ifstream ifile1("ThetaPolRcsPhi.out");
	}
	file1.close();
	}else{
		ofstream file("PhipolRcsTheta.out");
		for (int i = 0; i <= 900; i++) {
			file << RCSTheta[i] <<endl;
	//		file << RCSPhi[i] << endl;

			ifstream ifile("PhipolRcsTheta.out");
		}
		file.close();

		ofstream file1("PhipolRcsPhi.out");
		for (int i = 0; i <= 900; i++) {
	//		file << RCSTheta[i] <<endl;
			file1 << RCSPhi[i] << endl;
			ifstream ifile1("PhipolRcsPhi.out");
		}
		file1.close();

	}





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

//	delete[] Ec_antr;
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
	cout << (double)clock() / CLOCKS_PER_SEC << " S";
	return 0;

}

int GeoOpt(const Vector3 &raypoint, const BVH &bvh, Vector3 &hitpoint,
		Vector3 &k, Vector3 &E, float &t) {

	Vector3 phi_h, theta_h, Ei, Xc, Yc, Zc, theta_ci_h, phi_ci_h, theta_cr_h,
			phi_cr_h, Nz, Ed, m_h, Er;
	float theta_ci;
//	float phi_ci;
	int hitnumber = 0;
	Ei = E;
	t = 0;

	Ray ray(raypoint, k);
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

		if (theta_ci > 1e-6) {
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
			k = ray.d - 2 * (ray.d * nor) * nor;
			//    cout<<r.x<<","<<r.y<<","<<r.z<<"=k "<<endl;

			ray = Ray(I.hit, k);

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
			k = (-1) * ray.d;
			if (hitnumber % 2 == 1) {
				Er = (-1) * Ei;
			} else {
				Er = Ei;
			}
			ray = Ray(I.hit, k);
		}

		E = Er;
		hitpoint = I.hit;
		t += I.t;
//				    cout<<Ei.x<<","<<Ei.y<<","<<Ei.z<<" = ei"<<endl;
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
		const Vector3 &l, float &t) {
	// assuming vectors are all normalized
	float denom = n * l;
	if (denom > 1e-6) {
		Vector3 p0l0 = p0 - l0;
		t = p0l0 * n / denom;
		return (t >= 0);
	}

	return false;
}
void ComplexField(const Vector3 &R, const float &kr, const float frequency, CXV3 &C) {
	float multiply2 = 2.0;
	CXF expvalue = exp(j * (kr) / c * frequency * multiply2 * PI);
	C = CXV3(CXF(R.x, 0) * expvalue, CXF(R.y, 0) * expvalue,
			CXF(R.z, 0) * expvalue);
//	C = CXV3(CXF(R.x, 0) , CXF(R.y, 0) , CXF(R.z, 0) );

}

float Vector_Angle_d(const Vector3 &v1, const Vector3 &v2) {
	float l1 = length(v1);
	float l2 = length(v2);
	float angle = 0;
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

void rotate_THETA_then_PHI(Vector3 *R1, int n, float THETA, float PHI) {
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
void equation_twenty_seven(const CXV3 &E, const float &k0, const float cellSize, const Vector3 retpoint, const Vector3 &theta_h,
		const Vector3 &phi_h, CXF &Atheta, CXF &Aphi) {



//	float Sx = 0;
//	float Sy = 0;
//	float Xi = retpoint * phi_h;
//	float Yi = retpoint * theta_h;
//	cout<<Xi<<","<<Yi<<endl;
	float Ii = 1.0;
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
	if (Ex == CXF(0, 0) && Ey == CXF(0, 0)) {

		Atheta = CXF(0, 0);
		Aphi = CXF(0, 0);
	} else {

		Atheta = ( Ey ) * expvalue * cellSize * cellSize* (Ii);
//			cout<<"Atheta= "<<Atheta<<endl;
		Aphi = (Ex ) * expvalue * cellSize * cellSize * (Ii);
//			cout<<"Aphi= "<<Aphi<<endl;

	}

}

void GenrateAperture(const float&THETA, const float &PHI, const float &distance,
		const float cellSize, const int rnum,const int tnum,
		const float &Lx, const float &Ly, const int xdnum, const int ydnum, Vector3 *antr, Vector3 *antt,
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

//	for (int i = 0; i < rnum; i++) {
//		antr[i] = antr[i] - ant_Nd;
//	}

	for (int i = 0; i < tnum; i++) {

		antt[i] = antt[i] - ant_Nd;
	}


//	antc = antc - ant_Nd;

}

