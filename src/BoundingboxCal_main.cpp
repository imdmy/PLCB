#include "../include/BoundingBox.hpp"
// 只包含一些main()函数需要用到的头文件
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

int main(int argc, char **argv) {

  if (argc < 7) {
    std::cout << "Usage: " << argv[0] << " x1 y1 x2 y2 extendX extendY" << std::endl;
    std::cout << "你输入的点数不足" << std::endl;
    return 0;
  }
  std::vector<Point> OutputPts;
  std::vector<Point> points;

  double extendX = std::atof(argv[argc - 2]); // 首尾拓展距离
  double extendY = std::atof(argv[argc - 1]); // 侧边拓展距离

  std::cout << "extendX: " << extendX << " extendY: " << extendY << std::endl;

  // 调用时写成 vector输入
  for (int i = 1; i < argc - 2; i += 2) {
    Point p = {std::atof(argv[i]), std::atof(argv[i + 1])};
    points.push_back(p); // 将点存储起来
  }

  /*============两杆档距===============*/ // 可以将2塔和多塔融合到一个函数中，后续通过入参进行不同处理
  if (argc == 7) {                        // 两个点
    // 计算旋转包围框
    BoundingBox box = calculateRotatedBoundingBox(points[0], points[1]);
    // 计算扩展后的包围框顶点
    std::vector<Point> expandedCorners = calculateExpandedBoundingBox(
        box.centerPoint, box.width, box.height, box.angle, extendX, extendY);
    // 打印结果
    printCorners(expandedCorners);
    return 0;
    // 再根据包围框执行切割
    // do pdal
  }

  /*============多杆档距=============*/
  // 计算每条线段旋转角； 角度固定，旋转后可以按照中心偏移
  // 所有线段的偏移角； 获取的是线段与x轴的夹角。
  std::vector<double> rotate_Angle_all;
  for (int i{0}; i < points.size() - 1; i++) {
    // 获取基本旋转角 是否需要 (0.8*M_PI - s_Ag)定为旋转点
    double s_Ag{rotateAngle(points[i], points[i + 1])};
    rotate_Angle_all.push_back(s_Ag);
  }
  // 计算每三个点的夹角 计算的夹角总是基于 0 至 180°
  std::vector<double> angleBetweenPoints_all; // 保存所有线段组成的夹角
  for (int i = 0; i < points.size() - 2; i++) {
    double angle = angleBetweenPoints(points[i], points[i + 1],
                                      points[i + 2]); // 计算三点夹角
    angleBetweenPoints_all.push_back(angle);
  }

  std::vector<std::array<Point, 4>> waitRotatePoints; // 保存待旋转点
  // 计算每个线段的两个待旋转点进行绑定，可以绑定在待旋转点中，与原线段的两个点一起。
  for (int i{0}; i < points.size() - 2; i++) {
    std::array<Point, 4> rotatePt = widthBetweenAngle(
        points[i], points[i + 1], points[i + 2], angleBetweenPoints_all[i],
        extendY); // 单线段中由两个点生成两个的待旋转点。
    waitRotatePoints.push_back(rotatePt);
  }

  // 输出旋转点对，两个为一对，并进行排序
  std::vector<std::array<Point, 2>> RotatedPoints; // 保存旋转后的点
  // 后续四个点为组一起旋转
  // 分类判断 第一个点基于第二个点在第几象限  考虑如果位于象限上怎么办？
  // 如何处理？
  for (const auto &[index, unroPts] :
       std::ranges::views::enumerate(waitRotatePoints)) { // C23标准
    if ((unroPts[0].x - unroPts[1].x) >= 0 &&
        (unroPts[0].y - unroPts[1].y) >
            0) { // 第一个点位于基点第一象限 顺时针旋转
                 // 其实就是对夹角那两个点进行旋转，无需对基点（基准）坐标进行旋转
      bool direction{true};
      double rotateAngle = hfM_PI - rotate_Angle_all[index];
      std::array<Point, 2> rotatedPoints_Array =
          rotate4Points(unroPts, rotateAngle, direction);
      RotatedPoints.push_back(rotatedPoints_Array);
    }
    if ((unroPts[0].x - unroPts[1].x) > 0 &&
        (unroPts[0].y - unroPts[1].y) <=
            0) { // 第一个点位于基点第四象限 顺时针旋转
      bool direction{true};
      double rotateAngle = hfM_PI + rotate_Angle_all[index];
      cout << std::format("旋转角为: {}", rotateAngle * 180 / M_PI);
      std::array<Point, 2> rotatedPoints_Array =
          rotate4Points(unroPts, rotateAngle, direction);
      RotatedPoints.push_back(rotatedPoints_Array);
    }

    if ((unroPts[0].x - unroPts[1].x) < 0 &&
        (unroPts[0].y - unroPts[1].y) >=
            0) { // 第一个点位于基点第二象限 逆时针旋转
      bool direction{false};
      double rotateAngle = hfM_PI - rotate_Angle_all[index];
      std::array<Point, 2> rotatedPoints_Array =
          rotate4Points(unroPts, rotateAngle, direction);
      RotatedPoints.push_back(rotatedPoints_Array);
    }
    if ((unroPts[0].x - unroPts[1].x) <= 0 &&
        (unroPts[0].y - unroPts[1].y) <
            0) { // 第一个点位于基点第三象限 逆时针旋转
      bool direction{false};
      double rotateAngle = hfM_PI + rotate_Angle_all[index];
      std::array<Point, 2> rotatedPoints_Array =
          rotate4Points(unroPts, rotateAngle, direction);
      RotatedPoints.push_back(rotatedPoints_Array);
    }
  }

  cout << "\n\n旋转后的点为: ";
  for (const auto &[i, point] : std::ranges::views::enumerate(RotatedPoints)) {
    cout << std::format("\n第{}组点的坐标为：", i + 1);
    cout << std::format("\n\tPoint{}:({},{})", 1, point[0].x, point[0].y);
    cout << std::format("\n\tPoint{}:({},{})", 2, point[1].x, point[1].y);
  }

  // 排序，按照逆时针排序
  // 首先计算好首位点
  BoundingBox box1 = calculateRotatedBoundingBox(points[0], points[1]);
  BoundingBox box2 = calculateRotatedBoundingBox(
      points[points.size() - 2],
      points[points.size() -
             1]); // 存储3个数的索引值分别为0 1 2 所以为.size()-1
  std::vector<Point> expandedCorners1 = calculateExpandedBoundingBox(
      box1.centerPoint, box1.width, box1.height, box1.angle, extendX, extendY);
  // 打印结果
  std::cout << "\n第一个线段的包围框为:";
  printCorners(expandedCorners1);
  std::vector<Point> expandedCorners2 = calculateExpandedBoundingBox(
      box2.centerPoint, box2.width, box2.height, box2.angle, extendX, extendY);
  // 打印结果
  std::cout << "\n最后一个线段的包围框为:";
  printCorners(expandedCorners2);

  // expandedCorners1[0] 是切割框第一个点，expandedCorners1[3]
  // 是切割框最后一个点
  OutputPts.push_back(expandedCorners1[0]);

  for (const auto &[id, pt] : std::ranges::views::enumerate(RotatedPoints)) {
    OutputPts.push_back(RotatedPoints[id][0]);
  }
  OutputPts.push_back(expandedCorners2[1]);
  OutputPts.push_back(expandedCorners2[2]);
  // 倒序输入错误
  for (const auto &[id, pt] :
       std::ranges::views::enumerate(RotatedPoints | std::views::reverse)) {
    OutputPts.push_back(
        pt[1]); // std::views::reverse 反向指挥对pt进行反向，id索引还是正向的
  }
  OutputPts.push_back(expandedCorners1[3]); // 最后一个点

  cout << "\n\n包围框为: ";
  for (const auto &[i, point] : std::ranges::views::enumerate(OutputPts)) {
    cout << std::format("\n\tPoint{}: ({:.10f},{:.10f})", i, point.x, point.y);
  }

  return 0;
}
