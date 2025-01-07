#include "BoundingBox.hpp"
#include <array>
#include <cmath>
#include <iostream>

// ������ת����
void getRotationMatrix(double angle, double &cosTheta, double &sinTheta) {
  cosTheta = std::cos(angle);
  sinTheta = std::sin(angle);
}

// ��ת�� ��ʱ��anticlockwise
Point rotatePoint_anticlockwise(const Point &p, double cosTheta,
                                double &sinTheta) {
  return Point{p.x * cosTheta - p.y * sinTheta,
               p.x * sinTheta + p.y * cosTheta};
}

// ˳ʱ����ת
Point rotatePoint_clockwise(const Point &p, double cosTheta, double sinTheta) {
  return Point{p.x * cosTheta + p.y * sinTheta,
               -p.x * sinTheta + p.y * cosTheta};
}
// ����˳ʱ����ת����
void rotatePoint_clockwise(Point &p, double angle,
                           const Point &offset) { // offsetƫ����
  p.x = p.x * std::cos(angle) + p.y * std::sin(angle);
  p.y = -p.x * std::sin(angle) + p.y * std::cos(angle);
  p.x += offset.x;
  p.y += offset.y;
}
// ������ʱ����ת����
void rotatePoint_anticlockwise(Point &p, double angle, const Point &offset) {
  p.x = p.x * std::cos(angle) - p.y * std::sin(angle);
  p.y = p.x * std::sin(angle) + p.y * std::cos(angle);
  p.x += offset.x;
  p.y += offset.y;
}
// ���庯���� ����std::array<Point, 4> points �ĸ���
// �Լ���ת�ǣ�������ת���array�ĸ��� ���ȶ���������array rotatedPoints
// ���ȼ�����ת�ǵ�sin cos�� Ȼ�����ÿ���������ת��Ȼ����push �����µ����
std::array<Point, 2> rotate4Points(const std::array<Point, 4> &waitroPoints,
                                   const double &angle, const bool direction) {
  std::array<Point, 2> rotatedPoints{
      // ������ת��������нǵ�
      waitroPoints[2], waitroPoints[3] // ׼������ת��
  };
  if (direction) { // direction = true ˳ʱ��
    for (Point &point : rotatedPoints) {
      rotatePoint_clockwise(
          point, angle,
          waitroPoints[1]); // ������ δƫ�ƵĴ���ת�㣬��ת�ǣ� ƫ����
    }
  } else // direction = false ��ʱ��
  {
    for (Point &point : rotatedPoints) {
      rotatePoint_anticlockwise(point, angle, waitroPoints[1]);
    }
  }
  return rotatedPoints;
}

// ������չ��İ�Χ���ĸ�����
std::vector<Point> calculateExpandedBoundingBox(const Point &center,
                                                double width, double height,
                                                double angle, double x,
                                                double y) {
  std::vector<Point> corners(4);

  // ������չ��Ŀ�Ⱥ͸߶�
  double expandedWidth = width + 2 * x;
  double expandedHeight = height + 2 * y;

  // ������ת����
  double cosTheta, sinTheta;
  getRotationMatrix(angle, cosTheta, sinTheta);

  // ������ε��ĸ�����������ĵ��λ�� δ��תʱ�ĸ�����������ĵ��ƫ��
  std::vector<Point> cornersRel{{-expandedWidth / 2, -expandedHeight / 2},
                                {expandedWidth / 2, -expandedHeight / 2},
                                {expandedWidth / 2, expandedHeight / 2},
                                {-expandedWidth / 2, expandedHeight / 2}};

  // ��ת��ƽ�Ƶ����ĵ�
  for (size_t i = 0; i < corners.size(); ++i) {
    corners[i] = rotatePoint_anticlockwise(
        cornersRel[i], cosTheta, sinTheta); // todo: ˳ʱ����ʱ����Ҫ���
    corners[i].x += center.x;
    corners[i].y += center.y;
  }

  return corners;
}

// ��ӡ���ζ���
void printCorners(const std::vector<Point> &corners) {
  for (const auto &corner : corners) {
    std::cout << "\nPoint: (" << corner.x << ", " << corner.y << ")";
  }
}

// �����Եڶ�����Ϊ���㣬��һ����Ϊ��ת��
// �Ȳ�����˳ʱ�����ʱ�룬�Ծ���ֵ��ȡ�Ƕȣ�
// Ȼ���жϵ�һ����͵ڶ�������������ж���˳ʱ�뻹����ʱ��
double rotateAngle(const Point &p1, const Point &p2) {
  double rotateangle;
  // ���㷽�������ͽǶ�
  double dx = std::abs(p2.x - p1.x); // ��ȡ������
  double dy = std::abs(p2.y - p1.y); // ��ȡ������
  rotateangle = std::atan2(dy, dx);  // ��ȡ�����߶���x��н� ����
  cout << std::format("box.angle:{}��", rotateangle * (180 / M_PI));
  return rotateangle;
}
// ������ת��Χ��
BoundingBox calculateRotatedBoundingBox(const Point &p1, const Point &p2) {
  BoundingBox box;

  // �������ĵ�
  box.centerPoint.x = (p1.x + p2.x) / 2.0;
  box.centerPoint.y = (p1.y + p2.y) / 2.0;

  // // ���㷽�������ͽǶ�
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  box.angle = std::atan2(dy, dx); //todo: ����������β��Χ�����ת�Ƿ���-PI �� PI 
  cout << std::format("\nbox.angle:{}", box.angle); 

  /* box.angle = rotateAngle(p1, p2); */

  // �����������
  box.width = std::sqrt(dx * dx + dy * dy);
  box.height = 0; // ���

  return box;
}

// ��������н� ��
double angleBetweenPoints(const Point &A, const Point &B, const Point &C) {
  // ����A->B��C->B
  double vecBA_x = A.x - B.x;
  double vecBA_y = A.y - B.y;
  double vecBC_x = C.x - B.x;
  double vecBC_y = C.y - B.y;
  // ��������AB��CB�ĵ��
  double dotProduct = vecBA_x * vecBC_x + vecBA_y * vecBC_y;
  // ��������AB��CB�Ĳ��
  // double crossProduct = vecBA_x * vecBC_y - vecBA_y * vecBC_x;
  // BAx * BCy - BAy * BCx
  // ��������AB��CB�ĳ���
  double lengthAC = std::sqrt(vecBA_x * vecBA_x + vecBA_y * vecBA_y);
  double lengthCB = std::sqrt(vecBC_x * vecBC_x + vecBC_y * vecBC_y);
  // ����н�cos(��)
  double cosTheta = dotProduct / (lengthAC * lengthCB);
  // ʹ�÷����Һ���acos����нǦ�
  // ����нǷ�Χʼ�ջ����0-180��֮�䣬std::acos���ط�Χ�� 0 - PI
  double angle = std::acos(cosTheta);

  // ���Խ�����ת��Ϊ����ʾ
  // cout << std::format("�н�Ϊ: {}", angle = angle * 180.0 / M_PI);
  // ��δ�����ȶ�angle�����޸ģ��ڽ�������������return���ɶȡ�
  cout << std::format("\n�н�Ϊ: {}��", angle * 180.0 / M_PI);

  return angle;
}

// �������׼���������Ҫ������չ�ĶԵ��Ǵ���13���޻���24���ޣ�Ȼ�������洢�Ե��������ж��δ���
std::array<Point, 4> widthBetweenAngle(const Point &A, const Point &B,
                                       const Point &C, const double &angle,
                                       const double &extendY) {
  std::array<Point, 4> rotatePt{};
  Point pointBextend1, pointBextend2;

  // �������ж�˳ʱ�뻹����ʱ��;
  double vecBA_x = A.x - B.x;
  double vecBA_y = A.y - B.y;
  double vecBC_x = C.x - B.x;
  double vecBC_y = C.y - B.y;
  double crossProduct = vecBA_x * vecBC_y - vecBA_y * vecBC_x;

  // angle ����
  double halfAngle = angle * 0.5;
  // C++ ��׼��� std::sin �������ܻ����ƽǶȡ���ˣ���ʹ�� std::sin(angle)
  // ֮ǰ��Ҫȷ�� angle �ǻ����ơ�
  double distenceUp2{
      extendY /
      std::tan(
          halfAngle)}; // ������߽�����ԭ���ƫ�ƣ�Ӧ����ת���ټ���pointBd�����꣬�������ȼ�������ת��
  // ������������ ��pointB.y + distenceUp2; pointB.x + extendY�� �� ��pointB.y -
  // distenceUp2, pointB.x - extendY�� �ٽ�����ת����
  // �����˳ʱ�룬����13���ޣ��������ʱ�룬�Ǿ���24���ޡ�
  if (crossProduct > 0) { // ��ʱ�� ��չ��λ��24����
    pointBextend1 = {-extendY, distenceUp2};
    pointBextend2 = {extendY, -distenceUp2};
  } else { // ˳ʱ�� ��չ��λ��13����
    pointBextend1 = {-extendY, -distenceUp2};
    pointBextend2 = {extendY, distenceUp2};
  }

  rotatePt = {A, B, pointBextend1, pointBextend2};
  // �����ʾrotatePt��
  std::cout << "\nRotate Points:" << std::endl;
  for (std::size_t id{0}; id < rotatePt.max_size(); ++id) {
    std::cout << std::format("{}. ({},{})\n", id + 1, rotatePt[id].x,
                             rotatePt[id].y);
  }

  return rotatePt;
}
// ���Է�����ڽ�����⽻��
std::vector<Point> AnglePoint_raw(const double Angle, const double &cornerWidth,
                                  const Point &middlePoint) {
  double theta_rad{0.5 * Angle * M_PI / 180};
  Point resultinward, resultoutward;
  resultinward.x = middlePoint.x + cornerWidth * sin(theta_rad);
  resultinward.y = middlePoint.y + cornerWidth * cos(theta_rad);
  cout << "δ��ת���ڽ���Ϊ�� (" << resultinward.x << ", " << resultinward.y
       << ")";
  resultoutward.x = middlePoint.x - cornerWidth * sin(theta_rad);
  resultoutward.y = middlePoint.y - cornerWidth * cos(theta_rad);
  cout << "δ��ת���⽻�� Ϊ�� (" << resultoutward.x << ", " << resultoutward.y
       << ")" << endl;
  std::vector<Point> anglepoint{resultinward, resultoutward};
  double cosTheta, sinTheta;
  getRotationMatrix(Angle, cosTheta, sinTheta);
  for (Point &p : anglepoint) {
    p = rotatePoint_anticlockwise(p, cosTheta, sinTheta);
  }

  return anglepoint;
}