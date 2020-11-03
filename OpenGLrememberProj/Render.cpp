#include "Render.h"
#include <string>

#include <sstream>
#include <iostream>

#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>

#include "MyOGL.h"

#include "Camera.h"
#include "Light.h"
#include "Primitives.h"

#include "GUItextRectangle.h"

//---------------------------------
#include<vector>
using namespace std;
//---------------------------------

bool textureMode = true;
bool lightMode = true;
//---------------------------------
//��������� alpha ���������
int alpha = 0;
//����� ��������
bool texChange = true;
//----------------------------------

//����� ��� ��������� ������
class CustomCamera : public Camera
{
public:
	//��������� ������
	double camDist;
	//���� �������� ������
	double fi1, fi2;

	
	//������� ������ �� ���������
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//������� ������� ������, ������ �� ����� ��������, ���������� �������
	void SetUpCamera()
	{
		//�������� �� ������� ������ ������
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist*cos(fi2)*cos(fi1),
			camDist*cos(fi2)*sin(fi1),
			camDist*sin(fi2));

		if (cos(fi2) <= 0)
			normal.setCoords(0, 0, -1);
		else
			normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt()
	{
		//������� ��������� ������
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //������� ������ ������


//����� ��� ��������� �����
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//��������� ������� �����
		pos = Vector3(1, 1, 3);
	}

	
	//������ ����� � ����� ��� ���������� �����, ���������� �������
	void  DrawLightGhismo()
	{
		glDisable(GL_LIGHTING);

		
		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale*0.08;
		s.Show();
		
		if (OpenGL::isKeyPressed('G'))
		{
			glColor3d(0, 0, 0);
			//����� �� ��������� ����� �� ����������
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//������ ���������
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale*1.5;
			c.Show();
		}

	}

	void SetUpLight()
	{
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// ��������� ��������� �����
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// �������������� ����������� �����
		// ������� ��������� (���������� ����)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// ��������� ������������ �����
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// ��������� ���������� ������������ �����
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //������� �������� �����




//������ ���������� ����
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//������ ���� ������ ��� ������� ����� ������ ����
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//������� ���� �� ���������, � ����� ��� ����
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();

		double k = 0, x = 0, y = 0;
		if (r.direction.Z() == 0)
			k = 0;
		else
			k = (z - r.origin.Z()) / r.direction.Z();

		x = k*r.direction.X() + r.origin.X();
		y = k*r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02*dy);
	}

	
}

void mouseWheelEvent(OpenGL *ogl, int delta)
{

	if (delta < 0 && camera.camDist <= 1)
		return;
	if (delta > 0 && camera.camDist >= 100)
		return;

	camera.camDist += 0.01*delta;

}

void keyDownEvent(OpenGL *ogl, int key)
{
	if (key == 'L')
	{
		lightMode = !lightMode;
	}

	if (key == 'T')
	{
		textureMode = !textureMode;
	}

	if (key == 'R')
	{
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}

	if (key == 'F')
	{
		light.pos = camera.pos;
	}
	//-------------------------------------
	//��������� ����� ���������
	if (key == 'W') {
		alpha++;
		if (alpha == 3) {
			alpha = 0;
		}
	}
	//����� ��������
	if (key == 'Q') {
		texChange = !texChange;
	}
	//----------------------------------------
}

void keyUpEvent(OpenGL *ogl, int key)
{
	
}


//����� ��� ������ ������
vector<string> nametex_str = { "devil.bmp", "Rabi.bmp", "Home.bmp", "Home2.bmp" };
//����� ����������� ������� 
GLuint texId[4];

//����������� ����� ������ ��������
void initRender(OpenGL *ogl)
{
	//��������� �������

	//4 ����� �� �������� �������
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//��������� ������ ��������� �������
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//�������� ��������
	glEnable(GL_TEXTURE_2D);
	

	//������ ����������� ���������  (R G B)
	RGBTRIPLE *texarray;

	//������ ��������, (������*������*4      4, ���������   ����, �� ������� ������������ �� 4 ����� �� ������� �������� - R G B A)
	char *texCharArray;
	int texW, texH;
	for (int i = 0; i < nametex_str.size(); i++) {
		const char* nametex = nametex_str[i].c_str();
		OpenGL::LoadBMP(nametex, &texW, &texH, &texarray);
		OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

		//���������� �� ��� ��������
		glGenTextures(1, &texId[i]);
		//������ ��������, ��� ��� ����� ����������� � ���������, ����� ����������� �� ����� ��
		glBindTexture(GL_TEXTURE_2D, texId[i]);

		//��������� �������� � �����������, � ���������� ��� ������  ��� �� �����
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);

		//�������� ������
		free(texCharArray);
		free(texarray);

		//������� ����
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	}


	//������ � ���� ����������� � "������"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// ������������ �������� : �� ����� ����� ����� 1
	glEnable(GL_NORMALIZE);

	// ���������� ������������� ��� �����
	glEnable(GL_LINE_SMOOTH); 


	//   ������ ��������� ���������
	//  �������� GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  ������� � ���������� �������� ���������(�� ���������), 
	//                1 - ������� � ���������� �������������� ������� ��������       
	//                �������������� ������� � ���������� ��������� ����������.    
	//  �������� GL_LIGHT_MODEL_AMBIENT - ������ ������� ���������, 
	//                �� ��������� �� ���������
	// �� ��������� (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;
}

//---------------------------------------------------
class Point {
public:
	double x;
	double y;
	double z;
	Point(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	Point(double A[]) {
		this->x = A[0];
		this->y = A[1];
		this->z = A[2];
	}
	Point(const Point& A) {
		this->x = A.x;
		this->y = A.y;
		this->z = A.z;
	}
	void DrawPoint() {
		glVertex3d(this->x, this->y, this->z);
	}
	void Normal3d(int i = 1) {
		if (i > 0) {
			glNormal3d(x, y, z);
		}
		else {
			//������ ����������� �������
				x *= -1;
				y *= -1;
				z *= -1;
				glNormal3d(x, y, z);
		}
	}
	~Point() {

	}

};
class PointXY {
public:
	double x;
	double y;
	PointXY(double x, double y) {
		this->x = x;
		this->y = y;
	}
	void TexCoord2d() {
		glTexCoord2d(x, y);
	}
};
//������� ���������� ����� �������
static double SearchDistancePoints(Point A, Point B) {
	double X = B.x - A.x;
	double Y = B.y - A.y;
	double Z = B.z - A.z;
	X = pow(X, 2);
	Y = pow(Y, 2);
	Z = pow(Z, 2);
	double r = sqrt(X + Y + Z);
	return r;
}
//������������� ���������� ����� ������ � ������������ � �������� ��������
PointXY UpdatePoint(Point A, bool set = false, vector<Point> points = { Point(0, 0, 0) }, double r = 0) {
	static bool installation = false;
	static double width = 0;
	static double height = 0;
	static double min_x, min_y, max_x, max_y;
	static PointXY new_O(0, 0);
	PointXY newCoord(2, 2);
	if (set == true) {
		//���������� ������� ������ ��� ����� ����������
		min_x = points[0].x, min_y = points[0].y, max_x = points[0].x, max_y = points[0].y;
		for (int i = 0; i < points.size(); i++) {
			if (min_x > points[i].x) {
				min_x = points[i].x;
			}
			if (max_x < points[i].x) {
				max_x = points[i].x;
			}
			if (min_y > points[i].y) {
				min_y = points[i].y;
			}
			if (max_y < points[i].y) {
				max_y = points[i].y;
			}
		}

		//������ ����� � ������ �������������� ���� ������� ���� ������
		width = abs(min_x) + abs(max_x) + r;
		height = abs(min_y) + abs(max_y) + r;

		new_O.x = A.x;
		new_O.y = A.y;

		//��������� �����������
		installation = true;

		return newCoord;
	}
	if (installation) {
		//������ ����� ��������� ������������ 
		newCoord.x = (A.x - new_O.x) / width;
		newCoord.y = (A.y - new_O.y) / height;

		return newCoord;
	}
	newCoord.x = 9;
	newCoord.y = 9;
	return newCoord;
}
//����� ������ (A - ������ �������, B - ����� �������)
Point SearchVector(Point A, Point B) {
	return Point(B.x - A.x, B.y - A.y, B.z - A.z);
}
//��������� ������������
Point VectorProduct(Point vectorA, Point vectorB) {
	Point result(0, 0, 0);
	result.x = vectorA.y * vectorB.z - vectorB.y * vectorA.z;
	result.y = -1 * vectorA.x * vectorB.z + vectorB.x * vectorA.z;
	result.z = vectorA.x * vectorB.y - vectorB.x * vectorA.y;
	return result;
}
//����� ����� �������
static double SearchVectorLength(Point vector) {
	double length = sqrt(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2));
	return length;
}
//���� ������� � ������
Point SearchNormal(Point A, Point B, Point C, bool flag = true) {

	//����� �������� �� ��������������
	Point vector_b = SearchVector(B, C);
	Point vector_a = SearchVector(B, A);

	//��������� ������������
	Point NormalVector = VectorProduct(vector_a, vector_b);

	//����� ����� �������
	double length = SearchVectorLength(NormalVector);

	if (flag) {
		//������������ �������
		NormalVector.x = NormalVector.x / length;
		NormalVector.y = NormalVector.y / length;
		NormalVector.z = NormalVector.z / length;
	}
	return NormalVector;
}
//������� �������� �������
Point SearchMidpoint(Point A, Point B) {
	Point midpoint(0, 0, 0);
	midpoint.x = (A.x + B.x) / 2;
	midpoint.y = (A.y + B.y) / 2;
	midpoint.z = (A.z + B.z) / 2;
	return midpoint;
}
//����� ������ ����������. ���� �������� ��� ����� ������� �� ���
void FindCenterCircle(Point A, Point B, Point C, Point& Center, double& r) {
	Center.x = -0.5 * (A.y * (pow(B.x, 2) + pow(B.y, 2) - pow(C.x, 2) - pow(C.y, 2))
		+ B.y * (pow(C.x, 2) + pow(C.y, 2) - pow(A.x, 2) - pow(A.y, 2))
		+ C.y * (pow(A.x, 2) + pow(A.y, 2) - pow(B.x, 2) - pow(B.y, 2))) /
		(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));

	Center.y = 0.5 * (A.x * (pow(B.x, 2) + pow(B.y, 2) - pow(C.x, 2) - pow(C.y, 2))
		+ B.x * (pow(C.x, 2) + pow(C.y, 2) - pow(A.x, 2) - pow(A.y, 2))
		+ C.x * (pow(A.x, 2) + pow(A.y, 2) - pow(B.x, 2) - pow(B.y, 2))) /
		(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));

	r = SearchDistancePoints(Center, A);
}
//---------------------------------------------------
//����� ����� �� ���� ���������� (����������)
void SearchPointsCD(Point C, Point D, vector<Point>& pointsCD) {
	//���
	double step = 0.01;

	//����������
	//����� ���������� � ��� ������
	Point midpoint = SearchMidpoint(C, D);
	double r = SearchDistancePoints(C, midpoint);
	//������������ ����� ����������
	Point newPoint(0, 0, 0);
	for (double i = 0; i >= -180; i -= step) {
		newPoint.x = midpoint.x + r * cos(i * M_PI / 180);
		newPoint.y = midpoint.y + r * sin(i * M_PI / 180);
		pointsCD.push_back(newPoint);
	}
}
//����� ����� �� ���� ���������� (���������)
void SearchPointsAG(Point A, Point G, Point M, vector<Point>& pointsAG) {
	//���
	double step = 0.01;

	Point midpoint(0, 0, 0);
	double r = 0;
	FindCenterCircle(A, G, M, midpoint, r);
	//������������ ����� ���������
	Point newPoint(0, 0, 0);
	for (double i = -60.25; i <= 82.87; i += step) {
		newPoint.x = midpoint.x + r * cos(i * M_PI / 180);
		newPoint.y = midpoint.y + r * sin(i * M_PI / 180);
		pointsAG.push_back(newPoint);
	}
}
//������ �������
void sides(Point A, Point B, Point C, Point D, Point E, Point F, Point G, double z, vector<Point>& pointsCD, vector<Point>& pointsAG) {
	Point A1(0, 0, z), B1(-4, 0, z), C1(-1, -3, z), D1(4, -3, z), E1(3, 1, z), F1(4, 5, z), G1(-1, 5, z);
	glBegin(GL_QUADS);

	//�����
	SearchNormal(A, B, B1).Normal3d(-1);
	glTexCoord2d(0, 0);
	A.DrawPoint();
	glTexCoord2d(1, 0);
	B.DrawPoint();
	glTexCoord2d(1, 1);
	B1.DrawPoint();
	glTexCoord2d(0, 1);
	A1.DrawPoint();


	SearchNormal(B, C, C1).Normal3d(-1);
	glTexCoord2d(0, 0);
	B.DrawPoint();
	glTexCoord2d(1, 0);
	C.DrawPoint();
	glTexCoord2d(1, 1);
	C1.DrawPoint();
	glTexCoord2d(0, 1);
	B1.DrawPoint();

	SearchNormal(D, E, E1).Normal3d(-1);
	glTexCoord2d(0, 0);
	D.DrawPoint();
	glTexCoord2d(1, 0);
	E.DrawPoint();
	glTexCoord2d(1, 1);
	E1.DrawPoint();
	glTexCoord2d(0, 1);
	D1.DrawPoint();

	SearchNormal(E, F, F1).Normal3d(-1);
	glTexCoord2d(0, 0);
	E.DrawPoint();
	glTexCoord2d(1, 0);
	F.DrawPoint();
	glTexCoord2d(1, 1);
	F1.DrawPoint();
	glTexCoord2d(0, 1);
	E1.DrawPoint();

	SearchNormal(F, G, G1).Normal3d(-1);
	glTexCoord2d(0, 0);
	F.DrawPoint();
	glTexCoord2d(1, 0);
	G.DrawPoint();
	glTexCoord2d(1, 1);
	G1.DrawPoint();
	glTexCoord2d(0, 1);
	F1.DrawPoint();

	//������ ����������
	int size = pointsCD.size();
	int j = size;
	for (int i = 0; i < size - 1; i++) {
		
		Point A(pointsCD[i].x, pointsCD[i].y, 0);
		Point B(pointsCD[i + 1].x, pointsCD[i + 1].y, 0);
		Point A1(pointsCD[i].x, pointsCD[i].y, z);
		Point B1(pointsCD[i+1].x, pointsCD[i+1].y, z);

		double t1 = j / (size - 1.0);
		j = j - 1;
		double t2 = j / (size - 1.0);

		SearchNormal(A, B, B1).Normal3d();
		glTexCoord2d(t1, 0);
		A.DrawPoint();
		glTexCoord2d(t2, 0);
		B.DrawPoint();
		glTexCoord2d(t1, 1);
		B1.DrawPoint();
		glTexCoord2d(t2, 1);
		A1.DrawPoint();
	}

	//������ ���������
	size = pointsAG.size();
	j = size;
	for (int i = 0; i < size - 1; i++) {

		Point A(pointsAG[i].x, pointsAG[i].y, 0);
		Point B(pointsAG[i + 1].x, pointsAG[i + 1].y, 0);
		Point A1(pointsAG[i].x, pointsAG[i].y, z);
		Point B1(pointsAG[i + 1].x, pointsAG[i + 1].y, z);

		double t1 = j / (size - 1.0);
		j = j - 1;
		double t2 = j / (size - 1.0);

		SearchNormal(A, B, B1).Normal3d();
		glTexCoord2d(t1, 0);
		A.DrawPoint();
		glTexCoord2d(t2, 0);
		B.DrawPoint();
		glTexCoord2d(t1, 1);
		B1.DrawPoint();
		glTexCoord2d(t2, 1);
		A1.DrawPoint();
	}
	glEnd();
}

//������ ���
void bottom(Point A, Point B, Point C, Point D, Point E, Point F, Point G, double z, vector<Point> pointsCD, vector<Point> pointsAG) {
	glBegin(GL_TRIANGLES);

	SearchNormal(B, A, C).Normal3d(-1);
	UpdatePoint(B).TexCoord2d();
	B.DrawPoint();
	UpdatePoint(A).TexCoord2d();
	A.DrawPoint();
	UpdatePoint(C).TexCoord2d();
	C.DrawPoint();

	SearchNormal(D, A, E).Normal3d(-1);
	UpdatePoint(D).TexCoord2d();
	D.DrawPoint();
	UpdatePoint(A).TexCoord2d();
	A.DrawPoint();
	UpdatePoint(E).TexCoord2d();
	E.DrawPoint();
	
	int size = pointsCD.size();
	//����������� ����������
	for (int i = 0; i < size - 1; i++) {
		SearchNormal(pointsCD[i], A, pointsCD[i + 1]).Normal3d();
		UpdatePoint(pointsCD[i]).TexCoord2d();
		pointsCD[i].DrawPoint();
		UpdatePoint(A).TexCoord2d();
		A.DrawPoint();
		UpdatePoint(pointsCD[i + 1]).TexCoord2d();
		pointsCD[i+1].DrawPoint();
	}

	size = pointsAG.size();
	//����������� ���������
	for (int i = 0; i < size - 1; i++) {
		if (i > size / 2) {
			SearchNormal(F, pointsAG[i], pointsAG[i + 1]).Normal3d(-1);
			UpdatePoint(F).TexCoord2d();
			F.DrawPoint();
		}
		else {
			SearchNormal(E, pointsAG[i], pointsAG[i + 1]).Normal3d(-1);
			UpdatePoint(E).TexCoord2d();
			E.DrawPoint();
		}
		UpdatePoint(pointsAG[i]).TexCoord2d();
		pointsAG[i].DrawPoint();
		UpdatePoint(pointsAG[i + 1]).TexCoord2d();
		pointsAG[i + 1].DrawPoint();
	}

	SearchNormal(F, pointsAG[size / 2], E).Normal3d();
	UpdatePoint(F).TexCoord2d();
	F.DrawPoint();
	UpdatePoint(pointsAG[size / 2]).TexCoord2d();
	pointsAG[size / 2].DrawPoint();
	UpdatePoint(E).TexCoord2d();
	E.DrawPoint();

	glEnd();
}

//������ ����
void top(Point A, Point B, Point C, Point D, Point E, Point F, Point G, double z, vector<Point> pointsCD, vector<Point> pointsAG) {
	Point A1(0, 0, z), B1(-4, 0, z), C1(-1, -3, z), D1(4, -3, z), E1(3, 1, z), F1(4, 5, z), G1(-1, 5, z);
	
	//����������� ����� ���������
	switch (alpha) {
	case 1:
		//�������� ����� ����������
		glEnable(GL_BLEND);
		//������ ����� ��� ������������� ��������� � ���������
		glBlendFunc(GL_ONE, GL_ONE);
		break;
	case 2:
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glColor4d(0.7, 0.1, 0.1, 0.6);
		break;
	}
	glBegin(GL_TRIANGLES);

	SearchNormal(B1, A1, C1).Normal3d();
	UpdatePoint(B1).TexCoord2d();
	B1.DrawPoint();
	UpdatePoint(A1).TexCoord2d();
	A1.DrawPoint();
	UpdatePoint(C1).TexCoord2d();
	C1.DrawPoint();
	
	SearchNormal(D1, A1, E1).Normal3d();
	UpdatePoint(D1).TexCoord2d();
	D1.DrawPoint();
	UpdatePoint(A1).TexCoord2d();
	A1.DrawPoint();
	UpdatePoint(E1).TexCoord2d();
	E1.DrawPoint();

	int size = pointsCD.size();
	//����������� ����������
	pointsCD[0].z = z;
	for (int i = 0; i < size - 1; i++) {
		pointsCD[i + 1].z = z;
		SearchNormal(pointsCD[i], A1, pointsCD[i + 1]).Normal3d(-1);
		UpdatePoint(pointsCD[i]).TexCoord2d();
		pointsCD[i].DrawPoint();
		UpdatePoint(A1).TexCoord2d();
		A1.DrawPoint();
		UpdatePoint(pointsCD[i + 1]).TexCoord2d();
		pointsCD[i+1].DrawPoint();
	}

	size = pointsAG.size();
	//����������� ���������
	pointsAG[0].z = z;
	//����������� ���������
	for (int i = 0; i < size - 1; i++) {
		pointsAG[i+ 1].z = z;
		if (i > size / 2) {
			SearchNormal(F1, pointsAG[i], pointsAG[i + 1]).Normal3d();
			UpdatePoint(F1).TexCoord2d();
			F1.DrawPoint();
		}
		else {
			SearchNormal(E1, pointsAG[i], pointsAG[i + 1]).Normal3d();
			UpdatePoint(E1).TexCoord2d();
			E1.DrawPoint();
		}
		UpdatePoint(pointsAG[i]).TexCoord2d();
		pointsAG[i].DrawPoint();
		UpdatePoint(pointsAG[i + 1]).TexCoord2d();
		pointsAG[i + 1].DrawPoint();
	}

	SearchNormal(F1, pointsAG[size / 2], E1).Normal3d(-1);
	UpdatePoint(F1).TexCoord2d();
	F1.DrawPoint();
	UpdatePoint(pointsAG[size / 2]).TexCoord2d();
	pointsAG[size / 2].DrawPoint();
	UpdatePoint(E1).TexCoord2d();
	E1.DrawPoint();

	glEnd();

	//��������� ����������
	if (alpha != 0) {
		glDisable(GL_BLEND);
	}
}

void Render(OpenGL *ogl)
{

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);

	glEnable(GL_DEPTH_TEST);
	if (textureMode)
		glEnable(GL_TEXTURE_2D);

	if (lightMode)
		glEnable(GL_LIGHTING);

	//��������������
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//��������� ���������
	//GLfloat amb[] = { 0.2, 0.2, 0.1, 1. };
	GLfloat amb[] = { 112 / 255, 19 / 255, 61 / 255, 0.5 };
	GLfloat dif[] = { 0.4, 0.65, 0.5, 1. };
	GLfloat spec[] = { 0.9, 0.8, 0.3, 1. };
	GLfloat sh = 0.1f * 256;


	//�������
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	//��������
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	//����������
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec); \
		//������ �����
		glMaterialf(GL_FRONT, GL_SHININESS, sh);

	//���� ���� �������, ��� ����������� (����������� ���������)
	glShadeModel(GL_SMOOTH);
	//===================================
	//������� ���  

	//������ �����
	Point A(0, 0, 0), B(-4, 0, 0), C(-1, -3, 0), D(4, -3, 0), E(3, 1, 0), F(4, 5, 0), G(-1, 5, 0), M(1, 1, 0);
	double z = 1;

	vector<Point> pointsCD, pointsAG;
	//����� ����� ���� ����������
	SearchPointsCD(C, D, pointsCD);
	//����� ����� ���� ���������
	SearchPointsAG(A, G, M, pointsAG);

	//����� ��������
	if (texChange) {
		glBindTexture(GL_TEXTURE_2D, texId[2]);
	}
	else {
		glBindTexture(GL_TEXTURE_2D, texId[3]);
	}
	//������ ������� ������� 
	glColor3d(0.8, 0.1, 0.4);
	sides(A, B, C, D, E, F, G, z, pointsCD, pointsAG);

	vector<Point> points = { A, C, D, B, E, F, G };
	//���������� ��������� ��������� ��� ����������� ��������
	double r = SearchDistancePoints(SearchMidpoint(C, D), C);
	UpdatePoint(Point(-4, -5, 0), true, points, r);

	//����� ��������
	if (texChange) {
		glBindTexture(GL_TEXTURE_2D, texId[0]);
	}
	else {
		glBindTexture(GL_TEXTURE_2D, texId[1]);
	}
	//������ ���
	glColor3d(0.5, 0.2, 0.6);
	bottom(A, B, C, D, E, F, G, z, pointsCD, pointsAG);

	//����� ��������
	if (!texChange) {
		glBindTexture(GL_TEXTURE_2D, texId[0]);
	}
	else {
		glBindTexture(GL_TEXTURE_2D, texId[1]);
	}
	//������� ������
	glColor3d(1, 29 / 256, 0 / 256);
	top(A, B, C, D, E, F, G, z, pointsCD, pointsAG);

   //��������� ������ ������

	
	glMatrixMode(GL_PROJECTION);	//������ �������� ������� ��������. 
	                                //(���� ��������� ��������, ����� �� ������������.)
	glPushMatrix();   //��������� ������� ������� ������������� (������� ��������� ������������� ��������) � ���� 				    
	glLoadIdentity();	  //��������� ��������� �������
	glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 //������� ����� ������������� ��������

	glMatrixMode(GL_MODELVIEW);		//������������� �� �����-��� �������
	glPushMatrix();			  //��������� ������� ������� � ���� (��������� ������, ����������)
	glLoadIdentity();		  //���������� �� � ������

	glDisable(GL_LIGHTING);



	GuiTextRectangle rec;		   //������� ����� ��������� ��� ������� ������ � �������� ������.
	rec.setSize(300, 200);
	rec.setPosition(10, ogl->getHeight() - 200 - 10);


	std::stringstream ss;
	ss << "Q - ������� ��������" << std::endl;
	ss << "W - ������� ��� ������������" << std::endl;
	ss << "T - ���/���� �������" << std::endl;
	ss << "L - ���/���� ���������" << std::endl;
	ss << "F - ���� �� ������" << std::endl;
	ss << "G - ������� ���� �� �����������" << std::endl;
	ss << "G+��� ������� ���� �� ���������" << std::endl;
	ss << "�����. �����: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
	ss << "�����. ������: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
	ss << "��������� ������: R="  << camera.camDist << ", fi1=" << camera.fi1 << ", fi2=" << camera.fi2 << std::endl;
	
	rec.setText(ss.str().c_str());
	rec.Draw();

	glMatrixMode(GL_PROJECTION);	  //��������������� ������� �������� � �����-��� �������� �� �����.
	glPopMatrix();


	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

}