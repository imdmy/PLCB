#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cstdlib>  // For std::atof
#include <cassert>
#include <format>

#ifndef M_PI
const double M_PI = 3.14159265358979323846; // 定义 π 常量
const double hfM_PI = 0.5 * M_PI;
#endif
using std::cout,std::endl;

struct Point {
    double x, y;
};

struct BoundingBox {
    Point centerPoint;
    double width, height;
    double angle; // 旋转角度（弧度）
};

// 计算旋转矩阵
void getRotationMatrix(double angle, double& cosTheta, double& sinTheta) {
    cosTheta = std::cos(angle);
    sinTheta = std::sin(angle);
}

// 旋转点 逆时针anticlockwise
Point rotatePoint_anticlockwise(const Point& p, double cosTheta, double& sinTheta) {
    return Point{
        p.x * cosTheta - p.y * sinTheta,
        p.x * sinTheta + p.y * cosTheta
    };
}

//顺时针旋转
Point rotatePoint_clockwise(const Point& p, double cosTheta, double sinTheta) {
    return Point{
        p.x * cosTheta + p.y * sinTheta,
        -p.x * sinTheta + p.y * cosTheta
    };
}
// 重载顺时针旋转函数
void rotatePoint_clockwise( Point& p, double angle) {
    p.x = p.x * std::cos(angle) + p.y * std::sin(angle); 
    p.y = - p.x * std::sin(angle) + p.y * std::cos(angle);
}
// 定义函数： 输入std::array<Point, 4> points 四个点 以及旋转角，返回旋转后的array四个点
// 首先定义待输出的array rotatedPoints
// 首先计算旋转角的sin cos， 然后遍历每个点进行旋转，然后再push 进入新的输出
std::array <Point, 4> rotate4Points(const std::array<Point, 4>& waitroPoints, const double& angle) {
    std::array<Point, 4> rotatedPoints{ 
        waitroPoints[0],waitroPoints[1],waitroPoints[2],waitroPoints[3] //准备旋转点
    }; 
    for (Point point : rotatedPoints) {
        rotatePoint_clockwise(point, angle);
    }
}
// 计算扩展后的包围框四个顶点
std::vector<Point> calculateExpandedBoundingBox(const Point& center, double width, double height, double angle, double x, double y) {
    std::vector<Point> corners(4);

    // 计算扩展后的宽度和高度
    double expandedWidth = width + 2 * x;
    double expandedHeight = height + 2 * y;

    // 计算旋转矩阵
    double cosTheta, sinTheta;
    getRotationMatrix(angle, cosTheta, sinTheta);

    // 计算矩形的四个角相对于中心点的位置 未旋转时四个点相对于中心点的偏移
    std::vector<Point> cornersRel{
        { -expandedWidth / 2, -expandedHeight / 2 },
        { expandedWidth / 2, -expandedHeight / 2 },
        { expandedWidth / 2, expandedHeight / 2 },
        { -expandedWidth / 2, expandedHeight / 2 }
    };

    // 旋转并平移到中心点
    for (size_t i = 0; i < corners.size(); ++i) {
        corners[i] = rotatePoint_anticlockwise(cornersRel[i], cosTheta, sinTheta);
        corners[i].x += center.x;
        corners[i].y += center.y;
    }

    return corners;
}

// 打印矩形顶点
void printCorners(const std::vector<Point>& corners) {
    for (const auto& corner : corners) {
        std::cout << "Point: (" << corner.x << ", " << corner.y << ")\n";
    }
}

// 总是以第二个点为顶点，第一个点为旋转点
// 先不考虑顺时针和逆时针，以绝对值求取角度，
// 然后判断第一个点和第二个点的坐标来判断是顺时针还是逆时针
double rotateAngle(const Point& p1, const Point& p2) {
    double rotateangle;
    // 计算方向向量和角度
    double dx = std::abs(p2.x - p1.x); //获取横轴跨度
    double dy = std::abs(p2.y - p1.y); //获取纵轴跨度
    rotateangle = std::atan2(dy, dx); //获取的是线段与x轴夹角 弧度
    cout << std::format("box.angle：{}°", rotateangle * (180/M_PI));
    return rotateangle;
}
// 计算旋转包围框
BoundingBox calculateRotatedBoundingBox(const Point& p1, const Point& p2) {
    BoundingBox box;

    // 计算中心点
    box.centerPoint.x = (p1.x + p2.x) / 2.0;
    box.centerPoint.y = (p1.y + p2.y) / 2.0;

    // 计算方向向量和角度
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    box.angle = std::atan2(dy, dx);
    cout << std::format("box.angle：{}", box.angle);

    // 计算两点距离
    box.width = std::sqrt(dx * dx + dy * dy);
    box.height = 0; //宽度

    return box;
}

// 计算三点夹角 度
double angleBetweenPoints(const Point& A, const Point& B, const Point& C) {
    //计算A->B和C->B
    double vecAB_x = A.x - B.x;
    double vecAB_y = A.y - B.y;
    double vecCB_x = C.x - B.x;
    double vecCB_y = C.y - B.y;
    // 计算向量AB和CB的点积
    double dotProduct = vecAB_x * vecCB_x + vecAB_y * vecCB_y;
    // 计算向量AB和CB的长度
    double lengthAC = std::sqrt(vecAB_x * vecAB_x + vecAB_y * vecAB_y);
    double lengthCB = std::sqrt(vecCB_x * vecCB_x + vecCB_y * vecCB_y);
    // 计算夹角cos(θ)
    double cosTheta = dotProduct / (lengthAC * lengthCB);
    // 使用反余弦函数acos计算夹角θ
    // 计算夹角范围始终会介于0-180°之间，std::acos返回范围是 0 - PI
    double angle = std::acos(cosTheta);

    // 可以将弧度转换为度显示
    // cout << std::format("夹角为: {}", angle = angle * 180.0 / M_PI); 这段代码会先对angle进行修改，在进行输出，后面的return会变成度。
    cout << std::format("夹角为: {}", angle * 180.0 / M_PI);

    return angle;
}

std::array<Point, 4> widthBetweenAngle(const Point& pointA, const Point& pointB, const double& angle, const double& extendY) {
    std::array<Point, 4> rotatePt{};
    Point pointC, pointD;
    // angle 弧度
    double halfAngle = angle * 0.5;
    // C++ 标准库的 std::sin 函数接受弧度制角度。因此，在使用 std::sin(angle) 之前需要确保 angle 是弧度制。
    double distenceUp2{extendY/std::tan(halfAngle)}; 
    // 所以两个点是 （pointB.y + distenceUp2; pointB.x + extendY） 与 （pointB.y - distenceUp2, pointB.x - extendY） 再进行旋转即可
    pointC = {
        pointB.y + distenceUp2, pointB.x + extendY
    };
    pointD = {
        pointB.y - distenceUp2, pointB.x - extendY
    };
    rotatePt = {
        pointA,
        pointB,
        pointC,
        pointD    
    };
    //输出显示rotatePt：
    std::cout << "Rotate Points:" << std::endl;
    for (std::size_t id{ 0 }; id < rotatePt.max_size(); ++id) {
        std::cout << std::format("{}. ({},{})\n", id+1, rotatePt[id].x, rotatePt[id].y);
    }

    return rotatePt;
}
 // 可以分类成内交点和外交点 
std::vector<Point> AnglePoint_raw( const double Angle, const double& cornerWidth, const Point& middlePoint ) {
    double theta_rad{ 0.5*Angle * M_PI / 180 };
    Point resultinward, resultoutward;
    resultinward.x = middlePoint.x + cornerWidth * sin(theta_rad);
    resultinward.y = middlePoint.y + cornerWidth * cos(theta_rad);
    cout << "未旋转的内交点为： (" << resultinward.x << ", " << resultinward.y << ")";
    resultoutward.x = middlePoint.x - cornerWidth * sin(theta_rad);
    resultoutward.y = middlePoint.y - cornerWidth * cos(theta_rad);
    cout << "未旋转的外交点 为： (" << resultoutward.x << ", " << resultoutward.y << ")" << endl;
    std::vector<Point> anglepoint{ resultinward, resultoutward };
    double cosTheta, sinTheta;
    getRotationMatrix(Angle, cosTheta, sinTheta);
    for ( Point& p : anglepoint) {
        p = rotatePoint_anticlockwise(p, cosTheta, sinTheta);
    }

    return anglepoint;
}

/*
1. 首先将基准线设定为 基准旋转点 连接的上一个线段;
2. 将基准旋转点 至于原点（0，0）;
3. 计算基准线中心的旋转角；角度固定，旋转后可以按照中心偏移
4. 根据基准线和下一段线段的夹角开度（arctan > x or < x）确定夹角点的排序方式，首位对应排序；
5. 遍历到最后一个点后确定包围点，设定倒数两个点的包围框的后两个点，同理包围框开头是初始两个点的包围框上方两个点。;
6. 根据建立的点排序调取pdal切割算法，或者对pcpatch进行处理。;
*/
int main(int argc, char** argv) {
    // 解析输入点和扩展值
    std::vector<Point> points; // 保存输入点
    double extendX = std::atof(argv[argc - 2]); // 首尾拓展距离
    double extendY = std::atof(argv[argc - 1]); // 侧边拓展距离

    for (int i = 1; i < argc - 2; i += 2) {
        Point p = { std::atof(argv[i]), std::atof(argv[i + 1]) };
        points.push_back(p); //将点存储起来
    }
    /*============两个点===============*/
    if (argc == 7) { // 两个点
        // 计算旋转包围框
        BoundingBox box = calculateRotatedBoundingBox(points[0], points[1]);

        // 计算扩展后的包围框顶点
        std::vector<Point> expandedCorners = calculateExpandedBoundingBox(box.centerPoint, box.width, box.height, box.angle, extendX, extendY);

        // 打印结果
        printCorners(expandedCorners);
        return 0;
        // 再根据包围框执行切割
        // do pdal
    }

    /*============多点=============*/
    // 计算每条线段旋转角； 角度固定，旋转后可以按照中心偏移
    // 存储所有线段的偏移角； 获取的是线段与x轴的夹角。
    std::vector<double> rotate_Angle_all;
    for (int i{0}; i < points.size()-1; i += 1) {
        //获取基本旋转角 是否需要 (0.8*M_PI - s_Ag)定为旋转点
        double s_Ag{ rotateAngle(points[i], points[i+1])};
        rotate_Angle_all.push_back(s_Ag);
    }
    // 计算每三个点的夹角 计算的夹角总是基于 0 至 180°
    std::vector<double> angleBetweenPoints_all;
    for (int i = 0; i < points.size()-2; i += 1){
        double angle = angleBetweenPoints(points[i], points[i+1], points[i+2]);
        angleBetweenPoints_all.push_back(angle);
    }
    //#
    std::vector<std::array<Point, 4>> waitRotatePoints; //保存待旋转点
    // 计算每个线段的两个待旋转点进行绑定，可以绑定在待旋转点中，与原线段的两个点一起。
    for (int i{ 0 }; i < points.size() - 2; i += 1) {
        std::array<Point, 4> rotatePt = widthBetweenAngle(points[i], points[i+1], angleBetweenPoints_all[i], extendY); //单线段中由两个点生成两个的待旋转点。
        waitRotatePoints.push_back(rotatePt);
    }

    std::vector<std::array<Point, 4>> RotatedPoints; //保存旋转后的点
    // 后续四个点为组一起旋转
    // 分类判断 第一个点基于第二个点在第几象限
    for (const auto& unroPts : waitRotatePoints){
        //std::array<Point, 4> = 
        if ((unroPts[0].x - unroPts[1].x) > 0 && (unroPts[0].y - unroPts[1].y) > 0) {//第一个点位于基点第一象限 顺时针旋转
            double rotateAngle = hfM_PI - rotate_Angle_all[0];
            std::array<Point, 4> = //顺时针旋转四个点
        }

    }
    /*
    std::vector<Point> points1 = AnglePoint_raw(rotate_Angle_all[0], cornerWidth, points[1]); // 此处的角度应该是基准线段对于y轴的旋转角
    for (auto p : points1) {
        cout << std::format("旋转后的交点为：({},{})", p.x, p.y);
    }
    cout << std::sqrt(0.5 * 0.5 * 2);
    */

    return 0;
}