#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
//Библиотеки для OpenGL
#include <windows.h>
#include "GL/glut.h"
#include <gl\gl.h>
#include <gl\glu.h>

//Константы
#define PI 3.1415926535
#define EPS 0.01

//Количество шагов
#define STEPS_DRAW_CIRCLE 200

//Цвета
#define CHAIN_A_COLOR 0,0,0//Цвет исходной цепочки
#define CHAIN_B_COLOR 0,0,150//Цвет ответа
#define ALLOWED_AREA_COLOR 150,150,150
#define AREA_COLOR 0,0,0
#define AXE_COLOR 50,50,50
#define FINAL_POINT_COLOR 255,0,0
#define BG_COLOR 1,1,1

using namespace std;

//Начальная ширина и высота окна
GLint width = 1024, height =512;

//Класс Точка
class Point{
	public:
	//координаты
	double x;
	double y;
	конструкторы
	Point():x(0),y(0){}	
	Point(double xGet, double yGet):x(xGet),y(yGet){}
	//расстояние до другой точки
	double Distance(Point p) const {
		return sqrt((p.x-(*this).x)*(p.x-(*this).x)+(p.y-(*this).y)*(p.y-(*this).y));
	}
	//операторы
	const Point operator+(const Point& right) const {
        Point p(x+right.x,y+right.y);
		return p;
    }
	const Point operator-(const Point& right) const {
        Point p(x-right.x,y-right.y);
		return p;
    }
	const Point operator*(const double right) const {
        Point p(x*right,y*right);
		return p;
    }
	//рисуем точку
	void Draw(){
		glBegin(GL_POINTS);
			glVertex2d(x,y);
		glEnd();	
	}
};
//переопределяем операторы ввода и вывода
ostream& operator<<(ostream& out, const Point& print){
	out<<"("<<print.x<<","<<print.y<<")";
	return out;
}
istream& operator>>(istream& in, Point& read){
	in>>read.x>>read.y;	
	return in;
}

//Точка, к которую мы хотим перевести цепочку
Point P;

//Центры графиков
Point center1(width/6,height/2);
Point center2(width/2,height/2);
Point center3(5*width/6,height/2);

//Переводим радианы в градусы и градусы в радианы

double RadiansToGradus(double angleRadians);
double RadiansToGradus(double angleRadians)
{
	double angleGradus;
	angleGradus=angleRadians/PI*180;
	return angleGradus;
	
}

double GradusToRadians(double angleGradus);
double GradusToRadians(double angleGradus)
{
	double angleRadians;
	angleRadians=angleGradus/180*PI;
	return angleRadians;
	
}

//Рисует окружность заданного радиуса с центром в заданной точке
void DrawCircle(Point center, double radius);
void DrawCircle(Point center, double radius){
	
		GLint drawVertex[STEPS_DRAW_CIRCLE][2];
		double t;
		
		for(int i=0;i<STEPS_DRAW_CIRCLE;i++){
			t=2*PI/STEPS_DRAW_CIRCLE*i;
			drawVertex[i][0]=center.x+radius*cos(t);
			drawVertex[i][1]=center.y+radius*sin(t);
		}	

		glVertexPointer(2, GL_INT, 0, drawVertex);
		glEnableClientState(GL_VERTEX_ARRAY);
		glDrawArrays(GL_LINE_LOOP, 0, STEPS_DRAW_CIRCLE);
		glDisableClientState(GL_VERTEX_ARRAY);		
	
	
}

//Класс кольцо
class Ring{
	public:
	double smallRadius;
	double bigRadius;
	
	//конструкторы
	Ring():smallRadius(0),bigRadius(0){}
	Ring(double smallRadiusGet,double bigRadiusGet):smallRadius(smallRadiusGet),bigRadius(bigRadiusGet){}
	
	//рисуем кольцо	
	void Draw(Point center) const {
		
		DrawCircle(center,bigRadius);
		if(smallRadius!=0)
			DrawCircle(center,smallRadius);
		
	}
	
};

//Класс отрезок (или вектор)
class LineSegment{
	public:
	Point A;//начало
	Point B;//конец
	//конструкторы
	LineSegment(){
		Point p(0,0);
		A=p;
		B=p;		
	}	
	LineSegment(Point AGet, Point BGet):A(AGet),B(BGet){}
	//угол между векторами
	double Angle(LineSegment CD) const {
		Point C=CD.A;
		Point D=CD.B;
		return acos(((B.x-A.x)*(D.x-C.x)+(B.y-A.y)*(D.x-C.x))/(A.Distance(B)*C.Distance(D)));;
	}
	//рисуем отрезок
	void Draw(){
		glBegin(GL_LINES);
			glVertex2d(A.x,A.y);
			glVertex2d(B.x,B.y);
		glEnd();
		glBegin(GL_POINTS);
			glVertex2d(A.x,A.y);
			glVertex2d(B.x,B.y);
		glEnd();	
	}
};

//Вспомогательный класс, в котором передаём ответ
class AnswerPoints{
	public:
	Point A;
	Point B;
};

//Находим точки пересечения двух окружностей
AnswerPoints CirclesIntersection(Point C1, double r1, Point C2, double r2);
AnswerPoints CirclesIntersection(Point C1, double r1, Point C2, double r2)
{
	AnswerPoints answer;
	double d=C1.Distance(C2);
	if(d<=(r1+r2+0.01) && d>=r1+r2-0.01){
		answer.A=answer.B=C1+(C2-C1)*(r1/d);
	}else if(d>r1+r2 || (d<r1 && d+r2<r1) || (d<r2 && d+r1<r2) ){
		answer.A=answer.B=*(new Point(99999999,999999999));	
	}else{
		Point mainDirect=(C2-C1)*(1/d);
		Point normalDirect(mainDirect.y,-mainDirect.x);
		double otnoshenie=r1/r2*((r2*r2-d*d-r1*r1)/(2*r1*d))/((r1*r1-d*d-r2*r2)/(2*r2*d));
		Point intersectionLines=C1+mainDirect*d*(1/(1+1/otnoshenie));
		double cosb=(r1*r1-d*d-r2*r2)/(2*r2*d);
		double lengthIntersectionToAnswer=r2*sin(acos(cosb));

		answer.A=normalDirect*lengthIntersectionToAnswer+intersectionLines;
		answer.B=normalDirect*(-lengthIntersectionToAnswer)+intersectionLines;
	}
	
	return answer;
}

//Класс, с помощью которого перебираем все углы
class AnglesVector{
	public:
	int length;//длина вектора
	double * angles;//значения углов
	
	//конструкторы, деструктор и оператор =
	AnglesVector(): length(0),angles(NULL){}
	AnglesVector(int lengthGet,double * anglesGet):length(lengthGet){
		angles=new double[length];
		for(int i=0;i<length;i++){
			angles[i]=anglesGet[i];
		}
	}
	AnglesVector(int lengthGet):length(lengthGet){
		angles=new double[length];
		for(int i=0;i<length;i++){
			angles[i]=0;
		}
	}
	AnglesVector(AnglesVector const& copy):length(copy.length){
		angles=new double[length];
		for(int i=0;i<length;i++){
			angles[i]=copy.angles[i];
		}
	}
	AnglesVector& operator=(const AnglesVector& right) {
		length=right.length;
		if (angles!=NULL)delete[] angles;		
		angles=new double[length];
		for(int i=0;i<length;i++){
			angles[i]=right.angles[i];
		}
        return *this;
    }
	~AnglesVector(){
		delete[] angles;
	}
	
	//последний вектор, заканчиваем перебор
	bool IsEnd(){
		for(int i=0;i<length;i++)
			if(angles[i]<2*PI)
				return false;
		return true;
	}
	//следующий вектор
	void PlusStep(double step){
		if(IsEnd()==1)return;
		for(int i=length-1;i>=0;i--){
			if(angles[i]<2*PI){
				angles[i]+=step;//ШАГ
				return;
			}else{
				angles[i]=0;
			}
		}
	}

};
//определяем оператор вывода
ostream& operator<<(ostream& out, const AnglesVector& print){
	for(int i=0;i<print.length;i++)
		out<<print.angles[i]<<" ";
	return out;
}


//Класс цепочка
class Chain{
	public:
	int length;//длина цепочки
	Point * vertex;//вершины цепочки
	
	//конструкторы, деструкторы и оператор =
	Chain(): length(1){//(0,0)
		vertex=new Point[1];
	}
	Chain(int lengthGet,Point * vertexGet):length(lengthGet){
		vertex=new Point[length];
		for(int i=0;i<length;i++){
			vertex[i]=vertexGet[i];
		}
	}
	Chain(int lengthGet):length(lengthGet){
		vertex=new Point[length];
		Point p(0,0);
		for(int i=0;i<length;i++){
			vertex[i]=p;
		}
	}
	Chain(Chain const& copy):length(copy.length){
		vertex=new Point[length];
		for(int i=0;i<length;i++){
			vertex[i]=copy.vertex[i];
		}
	}
	Chain& operator=(const Chain& right) {
		length=right.length;
		if (vertex!=NULL)delete[] vertex;		
		vertex=new Point[length];
		for(int i=0;i<length;i++){
			vertex[i]=right.vertex[i];
		}
        return *this;
    }
	~Chain(){
		delete[] vertex;
	}
	
	//расстояние между i-ой и i+1-ой вершиной
	double DistanceVertexes(int i) const{
		if(!(i>=0 && (i+1)<length))
			return -9999999;
		return vertex[i].Distance(vertex[i+1]);		
	}
	//угол альфа_{i,i+1}
	double AngleVertexes(int i) const{
		if(!(i>=0 && (i+1)<length))
			return 99999999;
		if(i==0){
			Point O(0,0);
			Point A(50,0);
			LineSegment OA(O,A);
			LineSegment AB(O,vertex[1]);			
			return OA.Angle(AB);//Angle(O,B,O,vertex[i+1]);
		}else{
			LineSegment AB(vertex[i-1],vertex[i]);
			LineSegment BC(vertex[i],vertex[i+1]);			
			return AB.Angle(BC);//Angle(O,B,O,vertex[i+1]);
		}
		return 999999999;
	}
	//Область допустимых положений конца цепочки
	Ring AllowedArea(){
		Ring r;
		double sum=0;
		for(int i=1;i<length-1;i++){
			sum+=DistanceVertexes(i);
		}
		if(DistanceVertexes(0)-sum < 0)
			r.smallRadius=0;
		else
			r.smallRadius=DistanceVertexes(0)-sum;
		sum+=DistanceVertexes(0);
		r.bigRadius=sum;
		return r;
	}
	//Находится ли точка в области допустимый положений конца цепочки
	bool InAllowedArea(Point p){
		Ring r=AllowedArea();
		//cout<<"vertex[0].Distance(p)="<<vertex[0].Distance(p)<<endl;
		if(vertex[0].Distance(p)>=(r.smallRadius-EPS) && vertex[0].Distance(p)<=(r.bigRadius+EPS))
			return true;
		return false;
	}
	//Рисуем все области допустимых положений конца цепочки
	void DrawAllAllowedAreas(Point center){
		AllowedArea().Draw(center);
		for(int i=1;i<length;i++){
			Chain tmpChain(length-i,vertex);
			tmpChain.AllowedArea().Draw(center);
		}
		for(int i=length-1;i>=1;i--){
			Chain tmpChain(length-i,vertex+i-1);
			tmpChain.AllowedArea().Draw(center+vertex[i]);
		}
		
	}
	//Рисуем цепочку
	void Draw(Point center){
		LineSegment AB;
		for(int i=0;i<length-1;i++){
			AB.A=center+vertex[i];
			AB.B=center+vertex[i+1];
			AB.Draw();
		}
	}
	//Расстояние до цепочки
	double Distance(Chain toChain){
		double sum=0;
		for(int i=0;i<length;i++){
			sum+=vertex[i].Distance(toChain.vertex[i]);
		}
		return sum;
	}
	//Алгоритм А1
	Chain FindMinChainA1(double step){
		Chain minChain(length),tmpChain(length);
		if(!InAllowedArea(P))
			return minChain;
		AnglesVector angles(length-1-2);
		AnswerPoints intersectionCircles; 
		double minDistance=99999999;
		for(;!angles.IsEnd();angles.PlusStep(step)){
			tmpChain.vertex[length-1]=P;
			tmpChain.vertex[0]=vertex[0];
			for(int i=1;i<length-2;i++){
				tmpChain.vertex[i].x=tmpChain.vertex[i-1].x+DistanceVertexes(i-1)*cos(angles.angles[i-1]);
				tmpChain.vertex[i].y=tmpChain.vertex[i-1].y+DistanceVertexes(i-1)*sin(angles.angles[i-1]);
			}
			intersectionCircles=CirclesIntersection(tmpChain.vertex[length-3],DistanceVertexes(length-3),tmpChain.vertex[length-1],DistanceVertexes(length-2));
			if(intersectionCircles.A.Distance(vertex[length-2])<intersectionCircles.B.Distance(vertex[length-2]))
				tmpChain.vertex[length-2]=intersectionCircles.A;
			else
				tmpChain.vertex[length-2]=intersectionCircles.B;			
			if(Distance(tmpChain)<minDistance){
				minDistance=Distance(tmpChain);
				minChain=tmpChain;
			}			
		}		
		return minChain;
	}
	//Алгоритм А2
	Chain FindMinChainA2(){
		Chain minChain(length);
		if(!InAllowedArea(P))
			return minChain;
		minChain.vertex[length-1]=P;
		minChain.vertex[0]=vertex[0];
		if(length==2) return minChain;
		Ring allowedArea=AllowedArea();
		AnswerPoints allowedAreaIntersection,intersectionCircles;
		Point P1,P2,P3,P4,P5,napr,B,A;
		
		for(int i=length-2;i>=2;i--){
			Chain tmpChain(i+1,vertex);
			allowedArea=tmpChain.AllowedArea();
			allowedAreaIntersection=CirclesIntersection(vertex[0],allowedArea.bigRadius,minChain.vertex[i+1],DistanceVertexes(i));
			P1=allowedAreaIntersection.A;
			P2=allowedAreaIntersection.B;
			B=minChain.vertex[i+1];
			A=vertex[i+1];
			napr=(B-A)*(1/(A.Distance(B)));
			P5=B+napr*DistanceVertexes(i);
			P4=B-napr*DistanceVertexes(i);
			if(P4.Distance(vertex[i])<P5.Distance(vertex[i]))
				P3=P4;
			else
				P3=P5;			
			if(tmpChain.InAllowedArea(P3))//P3 в нужной зоне
				minChain.vertex[i]=P3;
			else{
				if(P1.Distance(vertex[i])<P2.Distance(vertex[i]))
					minChain.vertex[i]=P1;
				else
					minChain.vertex[i]=P2;
				
			}
		}
		//vertex 1
		intersectionCircles=CirclesIntersection(minChain.vertex[0],DistanceVertexes(0),minChain.vertex[2],DistanceVertexes(1));
			if(intersectionCircles.A.Distance(vertex[1])<intersectionCircles.B.Distance(vertex[1]))
				minChain.vertex[1]=intersectionCircles.A;
			else
				minChain.vertex[1]=intersectionCircles.B;			
		
		return minChain;

	}
	//Алгоритм А3	
	Chain FindMinChainA3(){
		Chain minChain(length);
		if(!InAllowedArea(P))
			return minChain;
		minChain.vertex[length-1]=P;
		minChain.vertex[0]=vertex[0];
		Ring r;
		int i,onPlace;
		for(i=length-2;i>=0;i--){
			Chain tmpChain(length-i,vertex+i);
			onPlace=i;
			r=tmpChain.AllowedArea();
			if(tmpChain.InAllowedArea(P))break;
		}
		//cout<<"onPlace="<<onPlace<<endl;
		for(i=0;i<=onPlace;i++){
			minChain.vertex[i]=vertex[i];			
		}
		Chain tmpChain(length-onPlace,vertex+onPlace);
		Chain tmpChainMin=tmpChain.FindMinChainA2();
		for(i=onPlace+1;i<length;i++){
			minChain.vertex[i]=tmpChainMin.vertex[i-onPlace];			
		}
		
	
		return minChain;

	}
	
		
};
//Переопределяем операторы ввода и вывода
ostream& operator<<(ostream& out, const Chain& print){
		
	out<<"Chain:"<<endl;
	out<<"	length="<<print.length<<" ("<<print.length-1<<" points)"<<endl;
	out<<"	";
	for(int i=0;i<print.length;i++)
		out<<print.vertex[i]<<" ";
	out<<endl;
	out<<"	d: ";
	for(int i=0;i<print.length-1;i++)
		out<<print.DistanceVertexes(i)<<" ";
	out<<endl;
	
	return out;
}
istream& operator>>(istream& in, Chain& read)
{
	in>>read.length;
	for(int i=0;i<read.length;i++)
		in>>read.vertex[i];
	return in;
}

//Определяем начальную цепочку и получаемые по алгоритмам
Chain firstChain;
Chain minChainA1;
Chain minChainA2;
Chain minChainA3;

//Перебираем все точки и строим область
void GoAllPoints(Point start,Point end,Point center);//start<end
void GoAllPoints(Point start,Point end,Point center){
	
	P=start;
	Point drawP;
	Point stepX(1,0),stepY(0,1);
		
	for(;P.x<end.x;P=P+stepX){
		for(;P.y<end.y;P=P+stepY){
			cout<<"P="<<P<<endl;
			if(!firstChain.InAllowedArea(P))continue;
			minChainA1=firstChain.FindMinChainA1(2*PI/500);//точнее!!
			minChainA3=firstChain.FindMinChainA3();
			if(minChainA1.Distance(minChainA3)<40){
				drawP=center+P;
				drawP.Draw();
			}
		}
		P.y=start.y;
		if(P.x>end.x)return;
	}
	
}

//Эта функция рисует всё
void Draw();
void Draw(){

    glClearColor(BG_COLOR, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
	//рисуем области допустимых положений конца цепочки
	glLineWidth(1);
	glColor3ub(ALLOWED_AREA_COLOR);
	firstChain.DrawAllAllowedAreas(center1);
	firstChain.DrawAllAllowedAreas(center2);
	firstChain.DrawAllAllowedAreas(center3);
	
	//рисуем оси
	glColor3ub(AXE_COLOR);
	glBegin(GL_LINES);
		glVertex2d(0,height/2);
		glVertex2d(width,height/2);

		glVertex2d(center1.x,0);
		glVertex2d(center1.x,height);

		glVertex2d(width/3,0);
		glVertex2d(width/3,height);

		glVertex2d(center2.x,0);
		glVertex2d(center2.x,height);

		glVertex2d(2*width/3,0);
		glVertex2d(2*width/3,height);

		glVertex2d(center3.x,0);
		glVertex2d(center3.x,height);
	glEnd();

	//рисуем цепочки
	glLineWidth(2);
	glPointSize(5);
	glColor3ub(CHAIN_A_COLOR);
	firstChain.Draw(center1);
	glColor3ub(CHAIN_B_COLOR);
	minChainA1.Draw(center1);
	minChainA2.Draw(center2);
	minChainA3.Draw(center3);

	//строим область
/*	
	glColor3ub(AREA_COLOR);
	glPointSize(2);
	Point start(-150,-150),end(150,150);
	GoAllPoints(start,end,center2);
*/

	//рисуем точку, в которой заканчиваются цепочки
	glColor3ub(FINAL_POINT_COLOR);
	Point drawP=center1+P;
	drawP.Draw();
	
	
	glFinish();
	
}


//эта функция управляет выводом на экран
void Display(void)
{

	Draw();

}

//Функция вызывается при изменении размеров окна
void Reshape(GLint w, GLint h)
{
    width = w;
    height = h;

    //устанавливаем размеры области отображения
    glViewport(0, 0, w, h);

    //ортографическая проекция
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, w, 0, h, -1.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

//функция, отвечающая за нажатие клавиш
void Keyboard(unsigned char key, int x, int y)
{
	#define ESCAPE '\033'

    if(key==ESCAPE)
        exit(0);
	
	return;
	
}

//нажатие кнопок мыши
void Mouse(int button, int state,int x, int y)
{
	Point clickedCoordinates(x,height-y);//Это надо исправить на настоящие координаты
	if(button == GLUT_LEFT_BUTTON && state==GLUT_DOWN){	
		if(clickedCoordinates.Distance(center1)<=firstChain.AllowedArea().bigRadius && clickedCoordinates.Distance(center1)>=firstChain.AllowedArea().smallRadius)
			P=clickedCoordinates-center1;
			//пересчитываем цепочки, находящиеся на минимальном расстоянии
			minChainA1=firstChain.FindMinChainA1(2*PI/500);
			minChainA2=firstChain.FindMinChainA2();
			minChainA3=firstChain.FindMinChainA3();
		
	}
	
	Draw();
}



int main(int argc, char *argv[]){

	//считываем точку и цепочку
	cin>>P;
	cin>>firstChain;
	//выводим цепочку
	cout<<firstChain<<endl;
	//считаем цепочки, находящиеся на минимальном расстоянии
	minChainA1=firstChain.FindMinChainA1(2*PI/500);
	cout<<minChainA1<<endl;
	minChainA2=firstChain.FindMinChainA2();
	cout<<minChainA2<<endl;
	minChainA3=firstChain.FindMinChainA3();
	cout<<minChainA3<<endl;
	
	//инициализируем графику
	glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB);
    glutInitWindowSize(width, height);
    glutCreateWindow("Chain");

    glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
    glutKeyboardFunc(Keyboard);
	glutMouseFunc(Mouse);
		
	glutMainLoop();

	return 0;
}

