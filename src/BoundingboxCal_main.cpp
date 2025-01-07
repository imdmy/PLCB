#include "../include/BoundingBox.hpp"
// ֻ����һЩmain()������Ҫ�õ���ͷ�ļ�
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

int main(int argc, char **argv) {

  if (argc < 7) {
    std::cout << "Usage: " << argv[0] << " x1 y1 x2 y2 extendX extendY" << std::endl;
    std::cout << "������ĵ�������" << std::endl;
    return 0;
  }
  std::vector<Point> OutputPts;
  std::vector<Point> points;

  double extendX = std::atof(argv[argc - 2]); // ��β��չ����
  double extendY = std::atof(argv[argc - 1]); // �����չ����

  std::cout << "extendX: " << extendX << " extendY: " << extendY << std::endl;

  // ����ʱд�� vector����
  for (int i = 1; i < argc - 2; i += 2) {
    Point p = {std::atof(argv[i]), std::atof(argv[i + 1])};
    points.push_back(p); // ����洢����
  }

  /*============���˵���===============*/ // ���Խ�2���Ͷ����ںϵ�һ�������У�����ͨ����ν��в�ͬ����
  if (argc == 7) {                        // ������
    // ������ת��Χ��
    BoundingBox box = calculateRotatedBoundingBox(points[0], points[1]);
    // ������չ��İ�Χ�򶥵�
    std::vector<Point> expandedCorners = calculateExpandedBoundingBox(
        box.centerPoint, box.width, box.height, box.angle, extendX, extendY);
    // ��ӡ���
    printCorners(expandedCorners);
    return 0;
    // �ٸ��ݰ�Χ��ִ���и�
    // do pdal
  }

  /*============��˵���=============*/
  // ����ÿ���߶���ת�ǣ� �Ƕȹ̶�����ת����԰�������ƫ��
  // �����߶ε�ƫ�ƽǣ� ��ȡ�����߶���x��ļнǡ�
  std::vector<double> rotate_Angle_all;
  for (int i{0}; i < points.size() - 1; i++) {
    // ��ȡ������ת�� �Ƿ���Ҫ (0.8*M_PI - s_Ag)��Ϊ��ת��
    double s_Ag{rotateAngle(points[i], points[i + 1])};
    rotate_Angle_all.push_back(s_Ag);
  }
  // ����ÿ������ļн� ����ļн����ǻ��� 0 �� 180��
  std::vector<double> angleBetweenPoints_all; // ���������߶���ɵļн�
  for (int i = 0; i < points.size() - 2; i++) {
    double angle = angleBetweenPoints(points[i], points[i + 1],
                                      points[i + 2]); // ��������н�
    angleBetweenPoints_all.push_back(angle);
  }

  std::vector<std::array<Point, 4>> waitRotatePoints; // �������ת��
  // ����ÿ���߶ε���������ת����а󶨣����԰��ڴ���ת���У���ԭ�߶ε�������һ��
  for (int i{0}; i < points.size() - 2; i++) {
    std::array<Point, 4> rotatePt = widthBetweenAngle(
        points[i], points[i + 1], points[i + 2], angleBetweenPoints_all[i],
        extendY); // ���߶��������������������Ĵ���ת�㡣
    waitRotatePoints.push_back(rotatePt);
  }

  // �����ת��ԣ�����Ϊһ�ԣ�����������
  std::vector<std::array<Point, 2>> RotatedPoints; // ������ת��ĵ�
  // �����ĸ���Ϊ��һ����ת
  // �����ж� ��һ������ڵڶ������ڵڼ�����  �������λ����������ô�죿
  // ��δ���
  for (const auto &[index, unroPts] :
       std::ranges::views::enumerate(waitRotatePoints)) { // C23��׼
    if ((unroPts[0].x - unroPts[1].x) >= 0 &&
        (unroPts[0].y - unroPts[1].y) >
            0) { // ��һ����λ�ڻ����һ���� ˳ʱ����ת
                 // ��ʵ���ǶԼн��������������ת������Ի��㣨��׼�����������ת
      bool direction{true};
      double rotateAngle = hfM_PI - rotate_Angle_all[index];
      std::array<Point, 2> rotatedPoints_Array =
          rotate4Points(unroPts, rotateAngle, direction);
      RotatedPoints.push_back(rotatedPoints_Array);
    }
    if ((unroPts[0].x - unroPts[1].x) > 0 &&
        (unroPts[0].y - unroPts[1].y) <=
            0) { // ��һ����λ�ڻ���������� ˳ʱ����ת
      bool direction{true};
      double rotateAngle = hfM_PI + rotate_Angle_all[index];
      cout << std::format("��ת��Ϊ: {}", rotateAngle * 180 / M_PI);
      std::array<Point, 2> rotatedPoints_Array =
          rotate4Points(unroPts, rotateAngle, direction);
      RotatedPoints.push_back(rotatedPoints_Array);
    }

    if ((unroPts[0].x - unroPts[1].x) < 0 &&
        (unroPts[0].y - unroPts[1].y) >=
            0) { // ��һ����λ�ڻ���ڶ����� ��ʱ����ת
      bool direction{false};
      double rotateAngle = hfM_PI - rotate_Angle_all[index];
      std::array<Point, 2> rotatedPoints_Array =
          rotate4Points(unroPts, rotateAngle, direction);
      RotatedPoints.push_back(rotatedPoints_Array);
    }
    if ((unroPts[0].x - unroPts[1].x) <= 0 &&
        (unroPts[0].y - unroPts[1].y) <
            0) { // ��һ����λ�ڻ���������� ��ʱ����ת
      bool direction{false};
      double rotateAngle = hfM_PI + rotate_Angle_all[index];
      std::array<Point, 2> rotatedPoints_Array =
          rotate4Points(unroPts, rotateAngle, direction);
      RotatedPoints.push_back(rotatedPoints_Array);
    }
  }

  cout << "\n\n��ת��ĵ�Ϊ: ";
  for (const auto &[i, point] : std::ranges::views::enumerate(RotatedPoints)) {
    cout << std::format("\n��{}��������Ϊ��", i + 1);
    cout << std::format("\n\tPoint{}:({},{})", 1, point[0].x, point[0].y);
    cout << std::format("\n\tPoint{}:({},{})", 2, point[1].x, point[1].y);
  }

  // ���򣬰�����ʱ������
  // ���ȼ������λ��
  BoundingBox box1 = calculateRotatedBoundingBox(points[0], points[1]);
  BoundingBox box2 = calculateRotatedBoundingBox(
      points[points.size() - 2],
      points[points.size() -
             1]); // �洢3����������ֵ�ֱ�Ϊ0 1 2 ����Ϊ.size()-1
  std::vector<Point> expandedCorners1 = calculateExpandedBoundingBox(
      box1.centerPoint, box1.width, box1.height, box1.angle, extendX, extendY);
  // ��ӡ���
  std::cout << "\n��һ���߶εİ�Χ��Ϊ:";
  printCorners(expandedCorners1);
  std::vector<Point> expandedCorners2 = calculateExpandedBoundingBox(
      box2.centerPoint, box2.width, box2.height, box2.angle, extendX, extendY);
  // ��ӡ���
  std::cout << "\n���һ���߶εİ�Χ��Ϊ:";
  printCorners(expandedCorners2);

  // expandedCorners1[0] ���и���һ���㣬expandedCorners1[3]
  // ���и�����һ����
  OutputPts.push_back(expandedCorners1[0]);

  for (const auto &[id, pt] : std::ranges::views::enumerate(RotatedPoints)) {
    OutputPts.push_back(RotatedPoints[id][0]);
  }
  OutputPts.push_back(expandedCorners2[1]);
  OutputPts.push_back(expandedCorners2[2]);
  // �����������
  for (const auto &[id, pt] :
       std::ranges::views::enumerate(RotatedPoints | std::views::reverse)) {
    OutputPts.push_back(
        pt[1]); // std::views::reverse ����ָ�Ӷ�pt���з���id�������������
  }
  OutputPts.push_back(expandedCorners1[3]); // ���һ����

  cout << "\n\n��Χ��Ϊ: ";
  for (const auto &[i, point] : std::ranges::views::enumerate(OutputPts)) {
    cout << std::format("\n\tPoint{}: ({:.10f},{:.10f})", i, point.x, point.y);
  }

  return 0;
}
