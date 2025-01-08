#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib> // For std::atof
#include <format>
#include <iostream>
#include <ranges>
#include <vector>

#ifndef M_PI
const double M_PI = 3.14159265358979323846; // 定义 π 常量
const double hfM_PI = 0.5 * M_PI;
#endif
using std::cout, std::endl;

struct Point {
  double x, y;
};

struct BoundingBox {
  Point centerPoint;
  double width, height;
  double angle; // 旋转角度（弧度）
};

// 计算旋转矩阵
void getRotationMatrix(double angle, double &cosTheta, double &sinTheta);

// 旋转点 逆时针anticlockwise
Point rotatePoint_anticlockwise(const Point &p, double cosTheta,
                                double &sinTheta);

// 顺时针旋转
Point rotatePoint_clockwise(const Point &p, double cosTheta, double sinTheta);

// 重载顺时针旋转函数
Point rotatePoint_clockwise(const Point &p, const double angle, const Point &offset);

// 重载逆时针旋转函数
Point rotatePoint_anticlockwise(const Point &p, const double angle, const Point &offset);

// 定义函数： 输入std::array<Point, 4> points 四个点
// 以及旋转角，返回旋转后的array四个点 首先定义待输出的array rotatedPoints
// 首先计算旋转角的sin cos， 然后遍历每个点进行旋转，然后再push 进入新的输出
std::array<Point, 2> rotate4Points(const std::array<Point, 4> &waitroPoints,
                                   const double &angle, const bool direction);

// 计算扩展后的包围框四个顶点
std::vector<Point> calculateExpandedBoundingBox(const Point &center,
                                                double width, double height,
                                                double angle, double x,
                                                double y);

// 打印矩形顶点
void printCorners(const std::vector<Point> &corners);

// 总是以第二个点为顶点，第一个点为旋转点
// 先不考虑顺时针和逆时针，以绝对值求取角度，
// 然后判断第一个点和第二个点的坐标来判断是顺时针还是逆时针
double rotateAngle(const Point &p1, const Point &p2);
// 计算旋转包围框
BoundingBox calculateRotatedBoundingBox(const Point &p1, const Point &p2);

// 计算三点夹角度
double angleBetweenPoints(const Point &A, const Point &B, const Point &C);

// 计算出标准距离后，仍需要考虑拓展的对点是处于13象限还是24象限，然后对这个存储对点容器进行二次处理。
std::array<Point, 4> widthBetweenAngle(const Point &A, const Point &B,
                                       const Point &C, const double &angle,
                                       const double &extendY);

// 可以分类成内交点和外交点
std::vector<Point> AnglePoint_raw(const double Angle, const double &cornerWidth,
                                  const Point &middlePoint);