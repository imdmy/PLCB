#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib> // For std::atof
#include <format>
#include <iostream>
#include <ranges>
#include <vector>

#ifndef M_PI
const double M_PI = 3.14159265358979323846; // ���� �� ����
const double hfM_PI = 0.5 * M_PI;
#endif
using std::cout, std::endl;

struct Point {
  double x, y;
};

struct BoundingBox {
  Point centerPoint;
  double width, height;
  double angle; // ��ת�Ƕȣ����ȣ�
};

// ������ת����
void getRotationMatrix(double angle, double &cosTheta, double &sinTheta);

// ��ת�� ��ʱ��anticlockwise
Point rotatePoint_anticlockwise(const Point &p, double cosTheta,
                                double &sinTheta);

// ˳ʱ����ת
Point rotatePoint_clockwise(const Point &p, double cosTheta, double sinTheta);

// ����˳ʱ����ת����
Point rotatePoint_clockwise(const Point &p, const double angle, const Point &offset);

// ������ʱ����ת����
Point rotatePoint_anticlockwise(const Point &p, const double angle, const Point &offset);

// ���庯���� ����std::array<Point, 4> points �ĸ���
// �Լ���ת�ǣ�������ת���array�ĸ��� ���ȶ���������array rotatedPoints
// ���ȼ�����ת�ǵ�sin cos�� Ȼ�����ÿ���������ת��Ȼ����push �����µ����
std::array<Point, 2> rotate4Points(const std::array<Point, 4> &waitroPoints,
                                   const double &angle, const bool direction);

// ������չ��İ�Χ���ĸ�����
std::vector<Point> calculateExpandedBoundingBox(const Point &center,
                                                double width, double height,
                                                double angle, double x,
                                                double y);

// ��ӡ���ζ���
void printCorners(const std::vector<Point> &corners);

// �����Եڶ�����Ϊ���㣬��һ����Ϊ��ת��
// �Ȳ�����˳ʱ�����ʱ�룬�Ծ���ֵ��ȡ�Ƕȣ�
// Ȼ���жϵ�һ����͵ڶ�������������ж���˳ʱ�뻹����ʱ��
double rotateAngle(const Point &p1, const Point &p2);
// ������ת��Χ��
BoundingBox calculateRotatedBoundingBox(const Point &p1, const Point &p2);

// ��������нǶ�
double angleBetweenPoints(const Point &A, const Point &B, const Point &C);

// �������׼���������Ҫ������չ�ĶԵ��Ǵ���13���޻���24���ޣ�Ȼ�������洢�Ե��������ж��δ���
std::array<Point, 4> widthBetweenAngle(const Point &A, const Point &B,
                                       const Point &C, const double &angle,
                                       const double &extendY);

// ���Է�����ڽ�����⽻��
std::vector<Point> AnglePoint_raw(const double Angle, const double &cornerWidth,
                                  const Point &middlePoint);