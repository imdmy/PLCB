#include "BoundingBox.hpp"
#include <array>
#include <cmath>
#include <iostream>

// 计算旋转矩阵
void getRotationMatrix(double angle, double &cosTheta, double &sinTheta) {
  cosTheta = std::cos(angle);
  sinTheta = std::sin(angle);
}

// 旋转点 逆时针anticlockwise
Point rotatePoint_anticlockwise(const Point &p, double cosTheta,
                                double &sinTheta) {
  return Point{p.x * cosTheta - p.y * sinTheta,
               p.x * sinTheta + p.y * cosTheta};
}

// 顺时针旋转
Point rotatePoint_clockwise(const Point &p, double cosTheta, double sinTheta) {
  return Point{p.x * cosTheta + p.y * sinTheta,
               -p.x * sinTheta + p.y * cosTheta};
}
// 重载顺时针旋转函数
Point rotatePoint_clockwise(const Point &p, const double angle, const Point &offset) { // offset偏移量
  Point returnP;
  returnP.x = p.x * std::cos(angle) + p.y * std::sin(angle);
  returnP.y = -p.x * std::sin(angle) + p.y * std::cos(angle);
  returnP.x += offset.x;
  returnP.y += offset.y;
  return returnP;
}

// 重载逆时针旋转函数
// 错误，第二个p.x已经使用了转换后的x坐标，导致旋转错误！！！
/* p.x = p.x * std::cos(angle) - p.y * std::sin(angle);
p.y = p.x * std::sin(angle) + p.y * std::cos(angle); */
Point rotatePoint_anticlockwise(const Point &p, const double angle, const Point &offset) {
  Point returnP;
  std::cout << "输入的偏移(即点B)为： " << std::format("{:.4f} {:.4f}", offset.x, offset.y)  << endl;
  std::cout << "旋转角" << angle * 180 / M_PI << "°" << endl;
  std::cout << " > 旋转的点为：" << std::format("{:.4f} {:.4f}", p.x, p.y) << endl;
  returnP.x = p.x * std::cos(angle) - p.y * std::sin(angle);
  returnP.y = p.x * std::sin(angle) + p.y * std::cos(angle);
  std::cout << "旋转后的坐标为(未偏移): " << std::format("{:.4f} {:.4f}", returnP.x, returnP.y) << endl;
  returnP.x += offset.x;
  returnP.y += offset.y;
  std::cout << "增加偏移后: " << std::format("{:.4f} {:.4f}", returnP.x, returnP.y) << endl;
  return returnP;
}
// 定义函数： 输入std::array<Point, 4> points 四个点
// 以及旋转角，返回旋转后的array四个点 首先定义待输出的array rotatedPoints
// 首先计算旋转角的sin cos， 然后遍历每个点进行旋转，然后再push 进入新的输出
// intput： {A, B, pointBextend1, pointBextend2}; 旋转角angle; 旋转方向direction;
// output: return 旋转后的内外夹点
std::array<Point, 2> rotate4Points(const std::array<Point, 4> &waitroPoints,
                                   const double &angle, const bool direction) {
  std::cout << "开始旋转, 目前旋转点" << waitroPoints[2].x << " "
            << waitroPoints[2].y << "; " << waitroPoints[3].x << " "
            << waitroPoints[3].y << endl;
            
/*   std::array<Point, 2> rotatedPoints{
      // 保存旋转后的两个夹角点
      waitroPoints[2], waitroPoints[3] // 准备待旋转点
  }; */

  std::array<Point, 2> rotatedPoints;
  if (direction == true) {  // direction = true 顺时针
    // 错误！std::array 不支持 push_back
    rotatedPoints[0] = rotatePoint_clockwise(waitroPoints[2], angle, waitroPoints[1]); // 分别对每个内外点进行旋转偏移
    rotatedPoints[1] = rotatePoint_clockwise(waitroPoints[3], angle, waitroPoints[1]);

  } else if (direction == false) // direction = false 逆时针
  {
    for (Point &point : rotatedPoints) {
      rotatedPoints[0] = rotatePoint_anticlockwise(waitroPoints[2], angle, waitroPoints[1]); // 参数： 未偏移的待旋转点，旋转角， 偏移量点B
      rotatedPoints[1] = rotatePoint_anticlockwise(waitroPoints[3], angle, waitroPoints[1]);
    }
  }
  std::cout << "旋转后的点坐标为： X: " << rotatedPoints[0].x << " "
            << rotatedPoints[0].y << "; Y: " << rotatedPoints[1].x << " "
            << rotatedPoints[1].y << endl;
  std::cin.get();
  return rotatedPoints;
}

// 计算扩展后的包围框四个顶点
std::vector<Point> calculateExpandedBoundingBox(const Point &center,
                                                double width, double height,
                                                double angle, double x,
                                                double y) {
  std::vector<Point> corners(4);

  // 计算扩展后的宽度和高度
  double expandedWidth = width + 2 * x;
  double expandedHeight = height + 2 * y;

  // 计算旋转矩阵
  double cosTheta, sinTheta;
  getRotationMatrix(angle, cosTheta, sinTheta);

  // 计算矩形的四个角相对于中心点的位置 未旋转时四个点相对于中心点的偏移
  std::vector<Point> cornersRel{{-expandedWidth / 2, -expandedHeight / 2},
                                {expandedWidth / 2, -expandedHeight / 2},
                                {expandedWidth / 2, expandedHeight / 2},
                                {-expandedWidth / 2, expandedHeight / 2}};

  // 旋转并平移到中心点
  for (size_t i = 0; i < corners.size(); ++i) {
    corners[i] = rotatePoint_anticlockwise(
        cornersRel[i], cosTheta, sinTheta); // todo: 顺时针逆时针需要解决
    corners[i].x += center.x;
    corners[i].y += center.y;
  }

  return corners;
}

// 打印矩形顶点
void printCorners(const std::vector<Point> &corners) {
  for (const auto &corner : corners) {
    std::cout << "\nPoint: (" << corner.x << ", " << corner.y << ")";
  }
}

// 总是以第二个点为顶点，第一个点为旋转点
// 先不考虑顺时针和逆时针，以绝对值求取角度，
// 然后判断第一个点和第二个点的坐标来判断是顺时针还是逆时针
double rotateAngle(const Point &p1, const Point &p2) {
  double rotateangle;
  // 计算方向向量和角度
  double dx = std::abs(p2.x - p1.x); // 获取横轴跨度
  double dy = std::abs(p2.y - p1.y); // 获取纵轴跨度
  rotateangle = std::atan2(dy, dx);  // 获取的是线段与x轴夹角 弧度
  cout << std::format("box.angle:{}°\t", rotateangle * (180 / M_PI));
  return rotateangle;
}

// 计算旋转包围框
BoundingBox calculateRotatedBoundingBox(const Point &p1, const Point &p2) {
  BoundingBox box;

  // 计算中心点
  box.centerPoint.x = (p1.x + p2.x) / 2.0;
  box.centerPoint.y = (p1.y + p2.y) / 2.0;

  // // 计算方向向量和角度
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  box.angle = std::atan2(dy, dx); //todo: 这里计算的首尾包围框的旋转角返回-PI 到 PI 
  // cout << std::format("\nbox.angle:{}", box.angle); 

  /* box.angle = rotateAngle(p1, p2); */

  // 计算两点距离
  box.width = std::sqrt(dx * dx + dy * dy);
  box.height = 0; // 宽度

  return box;
}

// 计算三点夹角 度
double angleBetweenPoints(const Point &A, const Point &B, const Point &C) {
  // 计算A->B和C->B
  double vecBA_x = A.x - B.x;
  double vecBA_y = A.y - B.y;
  double vecBC_x = C.x - B.x;
  double vecBC_y = C.y - B.y;
  // 计算向量AB和CB的点积
  double dotProduct = vecBA_x * vecBC_x + vecBA_y * vecBC_y;
  // 计算向量AB和CB的叉积
  // double crossProduct = vecBA_x * vecBC_y - vecBA_y * vecBC_x;
  // BAx * BCy - BAy * BCx
  // 计算向量AB和CB的长度
  double lengthAC = std::sqrt(vecBA_x * vecBA_x + vecBA_y * vecBA_y);
  double lengthCB = std::sqrt(vecBC_x * vecBC_x + vecBC_y * vecBC_y);
  // 计算夹角cos(θ)
  double cosTheta = dotProduct / (lengthAC * lengthCB);
  // 使用反余弦函数acos计算夹角θ
  // 计算夹角范围始终会介于0-180°之间，std::acos返回范围是 0 - PI
  double angle = std::acos(cosTheta);

  // 可以将弧度转换为度显示
  // cout << std::format("夹角为: {}", angle = angle * 180.0 / M_PI);
  // 这段代码会先对angle进行修改，在进行输出，后面的return会变成度。
  // cout << std::format("\n夹角为: {:.2f}°", angle * 180.0 / M_PI);

  return angle;
}

// 计算出标准距离后，仍需要考虑拓展的对点是处于13象限还是24象限，然后对这个存储对点容器进行二次处理。
std::array<Point, 4> widthBetweenAngle(const Point &A, const Point &B,
                                       const Point &C, const double &angle,
                                       const double &extendY) {
  std::array<Point, 4> rotatePt{};
  Point pointBextend1, pointBextend2;

  // 计算叉积判断顺时针还是逆时针;
  double vecBA_x = A.x - B.x;
  double vecBA_y = A.y - B.y;
  double vecBC_x = C.x - B.x;
  double vecBC_y = C.y - B.y;
  double crossProduct = vecBA_x * vecBC_y - vecBA_y * vecBC_x;

  // angle 弧度
  double halfAngle = angle * 0.5; // 中间点对是在角度的1/2处，所以需要除以2
  // C++ 标准库的 std::sin 函数接受弧度制角度。因此，在使用 std::sin(angle)
  // 之前需要确保 angle 是弧度制。
  double distenceUp2{
      extendY / std::tan(halfAngle) // 得出由于角度偏移计算出边界点基于原点的X轴距离偏移，应该旋转后，再加上pointBd的坐标，而不是先加上再旋转！
      }; 
  // 所以两个点是  外角（pointB.y + distenceUp2; pointB.x + extendY） 与 内角（pointB.y - distenceUp2, pointB.x - extendY） 再进行旋转即可
  // 如果是顺时针，就是13象限，如果是逆时针，那就是24象限。
  if (crossProduct > 0) { // 逆时针 点A在点B左侧 拓展点位于24象限
    pointBextend1 = {-extendY, distenceUp2};
    pointBextend2 = {extendY, -distenceUp2};
  } else { // 顺时针 点A在点B右侧 拓展点位于13象限
    pointBextend1 = {-extendY, -distenceUp2};
    pointBextend2 = {extendY, distenceUp2};
  }

  rotatePt = {A, B, pointBextend1, pointBextend2};
  // 输出显示rotatePt：
  std::cout << "\nRotate Points:" << std::endl;
  for (std::size_t id{0}; id < rotatePt.max_size(); ++id) {
    std::cout << std::format("{}. ({},{})\n", id + 1, rotatePt[id].x, rotatePt[id].y);
    std::cout << "内外夹点待旋转" << halfAngle * 180 / M_PI << "°\n" << std::endl;
  }

  return rotatePt;
}

// 可以分类成内交点和外交点
std::vector<Point> AnglePoint_raw(const double Angle, const double &cornerWidth,
                                  const Point &middlePoint) {
  double theta_rad{0.5 * Angle * M_PI / 180};
  Point resultinward, resultoutward;
  resultinward.x = middlePoint.x + cornerWidth * sin(theta_rad);
  resultinward.y = middlePoint.y + cornerWidth * cos(theta_rad);
  cout << "未旋转的内交点为： (" << resultinward.x << ", " << resultinward.y
       << ")";
  resultoutward.x = middlePoint.x - cornerWidth * sin(theta_rad);
  resultoutward.y = middlePoint.y - cornerWidth * cos(theta_rad);
  cout << "未旋转的外交点 为： (" << resultoutward.x << ", " << resultoutward.y
       << ")" << endl;
  std::vector<Point> anglepoint{resultinward, resultoutward};
  double cosTheta, sinTheta;
  getRotationMatrix(Angle, cosTheta, sinTheta);
  for (Point &p : anglepoint) {
    p = rotatePoint_anticlockwise(p, cosTheta, sinTheta);
  }

  return anglepoint;
}