#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>  // For std::atof
#include <cassert>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <format>

#ifndef M_PI
const double M_PI = 3.14159265358979323846; // 定义 π 常量
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
Point rotatePoint_anticlockwise(const Point& p, double cosTheta, double sinTheta) {
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
// 先不考虑顺时针和逆时针
double rotateAngle(const Point& p1, const Point& p2) {
    double rotateangle;
    // 计算方向向量和角度
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    rotateangle = std::atan2(dy, dx);
    cout << std::format("box.angle：{}", rotateangle);
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
    //计算AB和CB？
    double vecAB_x = A.x - B.x;
    double vecAB_y = A.y - B.y;
    double vecCB_x = C.x - B.x;
    double vecCB_y = C.y - B.y;
    /*
    // 计算向量A->C和B->C
    double vecAC_x = A.x - C.x;
    double vecAC_y = A.y - C.y;
    double vecBC_x = B.x - C.x;
    double vecBC_y = B.y - C.y;
    */
    
    // 计算向量AB和CB的点积
    double dotProduct = vecAB_x * vecCB_x + vecAB_y * vecCB_y;

    // 计算向量AB和CB的长度
    double lengthAC = std::sqrt(vecAB_x * vecAB_x + vecAB_y * vecAB_y);
    double lengthBC = std::sqrt(vecCB_x * vecCB_x + vecCB_y * vecCB_y);

    // 计算夹角cos(θ)
    double cosTheta = dotProduct / (lengthAC * lengthBC);

    // 使用反余弦函数acos计算夹角θ
    double angle = std::acos(cosTheta);

    // 将弧度转换为度
    angle = angle * 180.0 / M_PI;

    return angle;
}

double widthBetweenAngle(const Point& pointB, double angle, double width) {
    double halfAngle = angle * 0.5;
    std::cout << "对称角度: " << halfAngle;
    // C++ 标准库的 std::sin 函数接受弧度制角度。因此，在使用 std::sin(angle) 之前需要确保 angle 是弧度制。
    double angleInDegrees{ halfAngle };
    double angleInRadians = angleInDegrees * ( M_PI / 180);
    double cornerWidth = width / std::sin(angleInRadians); //sin() 弧度制下的对称角
    assert(cornerWidth >= width);
    std::cout << "旋转角扩展距离: " << cornerWidth << std::endl;
    std::cout << "直角邻边: " << cornerWidth << std::endl;
    return cornerWidth;
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
3. 计算基准线中心的旋转角；
4. 根据基准线和下一段线段的夹角开度（arctan > x or < x）确定夹角点的排序方式，首位对应排序；
5. 遍历到最后一个点后确定包围点，设定倒数两个点的包围框的后两个点，同理包围框开头是初始两个点的包围框上方两个点。;
6. 根据建立的点排序调取pdal切割算法，或者对pcpatch进行处理。;
*/
int main(int argc, char** argv) {
    // 解析输入点和扩展值
    std::vector<Point> points;
    double extendX = std::atof(argv[argc - 2]); // 首尾拓展距离
    double extendY = std::atof(argv[argc - 1]); // 侧边拓展距离

    for (int i = 1; i < argc - 2; i += 2) {
        Point p = { std::atof(argv[i]), std::atof(argv[i + 1]) };
        points.push_back(p);
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
    //存储所有线段的偏移角；
    std::vector<double> rotate_Angle;
    for (int i{0}; i < points.size() -2 ; i += 2) {
        //获取基本旋转角
        double s_Ag{ rotateAngle(points[i], points[i+1])};
        rotate_Angle.push_back(s_Ag);
    }

    double angle = angleBetweenPoints(points[0], points[1], points[2]);
    std::cout << angle << std::endl;

    double cornerWidth = widthBetweenAngle(points[1], angle, extendY);
    std::vector<Point> points1 = AnglePoint_raw(rotate_Angle[0], cornerWidth, points[1]); // 此处的角度应该是基准线段对于y轴的旋转角
    for (auto p : points1) {
        cout << std::format("旋转后的交点为：({},{})", p.x, p.y);
    }
    cout << std::sqrt(0.5 * 0.5 * 2);
    return 0;
}