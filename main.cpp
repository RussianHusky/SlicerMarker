#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <fstream>

/**
 * @brief 2D point.
 */
struct Point_2 {
    double x{};
    double y{};
};

/**
 * @brief Line segment defined by two points.
 */
struct Line_2 {
    Point_2 p1;
    Point_2 p2;
};

/**
 * @brief Axis-aligned bounding box.
 */
struct BoundingBox {
    double xmin{};
    double xmax{};
    double ymin{};
    double ymax{};
};

/**
 * @brief Compute bounding box of a set of points.
 */
BoundingBox computeBoundingBox(const std::vector<Point_2>& pts) {
    if (pts.empty()) {
        return {0.0, 0.0, 0.0, 0.0};
    }

    double xmin = pts[0].x;
    double xmax = pts[0].x;
    double ymin = pts[0].y;
    double ymax = pts[0].y;

    for (const auto& p : pts) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }

    return {xmin, xmax, ymin, ymax};
}

/**
 * @brief Remove duplicate points from a list using a distance tolerance.
 */
void uniquePoints(std::vector<Point_2>& pts, double eps = 1e-9) {
    std::vector<Point_2> result;
    for (const auto& p : pts) {
        bool dup = false;
        for (const auto& q : result) {
            double dx = p.x - q.x;
            double dy = p.y - q.y;
            if (dx * dx + dy * dy <= eps * eps) {
                dup = true;
                break;
            }
        }
        if (!dup) {
            result.push_back(p);
        }
    }
    pts = std::move(result);
}

/**
 * @brief Compute intersection of a line (n·p = c) with axis-aligned rectangle.
 *
 * Line is defined by normal vector n and scalar c:
 *      n.x * x + n.y * y = c
 *
 * @param n      Normal vector (unit length is preferred).
 * @param c      Constant term for the line equation.
 * @param bbox   Bounding box of the rectangle.
 * @param cosA   cos(angle) of hatch direction.
 * @param sinA   sin(angle) of hatch direction.
 * @return vector of intersection points (0, 1, 2 or more if passing through corners).
 */
std::vector<Point_2> clipLineWithRectangle(const Point_2& n,
                                           double c,
                                           const BoundingBox& bbox,
                                           double cosA,
                                           double sinA)
{
    const double eps = 1e-9;
    std::vector<Point_2> intersections;

    Point_2 d{cosA, sinA};

    if (std::fabs(n.y) > eps) {
        {
            double x = bbox.xmin;
            double y = (c - n.x * x) / n.y;
            if (y >= bbox.ymin - eps && y <= bbox.ymax + eps) {
                intersections.push_back({x, y});
            }
        }
        {
            double x = bbox.xmax;
            double y = (c - n.x * x) / n.y;
            if (y >= bbox.ymin - eps && y <= bbox.ymax + eps) {
                intersections.push_back({x, y});
            }
        }
    }

    if (std::fabs(n.x) > eps) {
        {
            double y = bbox.ymin;
            double x = (c - n.y * y) / n.x;
            if (x >= bbox.xmin - eps && x <= bbox.xmax + eps) {
                intersections.push_back({x, y});
            }
        }
        {
            double y = bbox.ymax;
            double x = (c - n.y * y) / n.x;
            if (x >= bbox.xmin - eps && x <= bbox.xmax + eps) {
                intersections.push_back({x, y});
            }
        }
    }

    uniquePoints(intersections);

    if (intersections.size() <= 2) {
        return intersections;
    }

    std::sort(intersections.begin(), intersections.end(),
              [&d](const Point_2& a, const Point_2& b) {
                  double ta = d.x * a.x + d.y * a.y;
                  double tb = d.x * b.x + d.y * b.y;
                  return ta < tb;
              });

    std::vector<Point_2> trimmed;
    trimmed.push_back(intersections.front());
    trimmed.push_back(intersections.back());
    return trimmed;
}

/**
 * @brief Generate hatch lines inside a rectangular contour.
 *
 * Assumptions:
 *  - contour describes axis-aligned rectangle (or square),
 *  - lines are infinite, then clipped to the rectangle bounds.
 *
 * @param contour Rectangular contour points.
 * @param angle_deg Angle of hatch (degrees, 0..180).
 * @param step Distance between parallel lines (same units as coordinates).
 * @return vector of hatch line segments.
 */
std::vector<Line_2> generateHatch(const std::vector<Point_2>& contour,
                                  double angle_deg,
                                  double step)
{
    std::vector<Line_2> result;

    if (contour.empty() || step <= 0.0) {
        return result;
    }

    const double pi = std::acos(-1.0);
    double angle_rad = angle_deg * pi / 180.0;

    double cosA = std::cos(angle_rad);
    double sinA = std::sin(angle_rad);

    Point_2 n{-sinA, cosA};

    BoundingBox bbox = computeBoundingBox(contour);

    std::vector<Point_2> corners = {
        {bbox.xmin, bbox.ymin},
        {bbox.xmax, bbox.ymin},
        {bbox.xmax, bbox.ymax},
        {bbox.xmin, bbox.ymax}
    };

    double min_c = std::numeric_limits<double>::infinity();
    double max_c = -std::numeric_limits<double>::infinity();

    for (const auto& p : corners) {
        double val = n.x * p.x + n.y * p.y;
        min_c = std::min(min_c, val);
        max_c = std::max(max_c, val);
    }

    if (!std::isfinite(min_c) || !std::isfinite(max_c)) {
        return result;
    }

    const double eps = 1e-9;

    for (double c = min_c; c <= max_c + eps; c += step) {
        auto pts = clipLineWithRectangle(n, c, bbox, cosA, sinA);
        if (pts.size() == 2) {
            result.push_back({pts[0], pts[1]});
        }
    }

    return result;
}

/**
 * @brief Export contour and hatch lines to a simple SVG file.
 *
 * @param contour Rectangular contour.
 * @param lines Hatch lines.
 * @param filename Path to SVG file.
 */
void exportToSvg(const std::vector<Point_2>& contour,
                 const std::vector<Line_2>& lines,
                 const std::string& filename)
{
    if (contour.empty()) {
        return;
    }

    BoundingBox bbox = computeBoundingBox(contour);

    double scale = 20.0;
    double padding = 20.0;

    double width  = (bbox.xmax - bbox.xmin) * scale + 2.0 * padding;
    double height = (bbox.ymax - bbox.ymin) * scale + 2.0 * padding;

    auto toScreenX = [=](double x) {
        return padding + (x - bbox.xmin) * scale;
    };

    auto toScreenY = [=](double y) {
        return padding + (bbox.ymax - y) * scale;
    };

    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Failed to open SVG file: " << filename << "\n";
        return;
    }

    ofs << R"(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)" << "\n";
    ofs << R"(<svg xmlns="http://www.w3.org/2000/svg" version="1.1")";
    ofs << " width=\"" << width << "\" height=\"" << height << "\"";
    ofs << " viewBox=\"0 0 " << width << " " << height << "\">" << "\n";

    ofs << R"(<rect x="0" y="0" width=")" << width << R"(" height=")" << height
        << R"(" fill="white" stroke="none"/>)" << "\n";

    ofs << R"(<polygon fill="none" stroke="black" stroke-width="1" points=")";
    for (const auto& p : contour) {
        ofs << toScreenX(p.x) << "," << toScreenY(p.y) << " ";
    }
    ofs << "\"/>\n";

    // Draw hatch lines
    ofs << R"(<g stroke="blue" stroke-width="0.8">)" << "\n";
    for (const auto& l : lines) {
        ofs << R"(<line x1=")" << toScreenX(l.p1.x)
            << R"(" y1=")" << toScreenY(l.p1.y)
            << R"(" x2=")" << toScreenX(l.p2.x)
            << R"(" y2=")" << toScreenY(l.p2.y)
            << R"(" />)" << "\n";
    }
    ofs << "</g>\n";

    ofs << "</svg>\n";
}

/**
 * @brief Simple command line arguments parser for --angle, --step, --svg.
 */
struct CmdOptions {
    double angle_deg = 45.0;
    double step      = 1.0;
    std::optional<std::string> svgFile;
};

CmdOptions parseCommandLine(int argc, char** argv) {
    CmdOptions opt;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg(argv[i]);

        if (arg == "--angle" && i + 1 < argc) {
            opt.angle_deg = std::stod(argv[++i]);
        } else if (arg == "--step" && i + 1 < argc) {
            opt.step = std::stod(argv[++i]);
        } else if (arg == "--svg" && i + 1 < argc) {
            opt.svgFile = std::string(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0]
                      << " --angle <deg> --step <value> [--svg <file>] < input\n";
            std::exit(0);
        }
    }

    return opt;
}

/**
 * @brief Read contour points from stdin.
 *
 * Format:
 *  N
 *  x0 y0
 *  x1 y1
 *  ...
 */
std::vector<Point_2> readContourFromStdin() {
    std::vector<Point_2> contour;

    std::size_t n = 0;
    if (!(std::cin >> n)) {
        return contour;
    }

    contour.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        Point_2 p{};
        if (!(std::cin >> p.x >> p.y)) {
            contour.clear();
            break;
        }
        contour.push_back(p);
    }

    return contour;
}

/**
 * @brief Entry point.
 */
int main(int argc, char** argv) {
    CmdOptions options = parseCommandLine(argc, argv);

    std::vector<Point_2> contour = readContourFromStdin();
    if (contour.empty()) {
        contour = {
            {0.0, 0.0},
            {10.0, 0.0},
            {10.0, 10.0},
            {0.0, 10.0}
        };
        std::cerr << "No valid contour read from stdin. "
                     "Using default square [0,10] x [0,10].\n";
    }

    auto lines = generateHatch(contour, options.angle_deg, options.step);

    std::cout << std::fixed << std::setprecision(3);
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const auto& l = lines[i];
        std::cout << "Line " << (i + 1) << ": ("
                  << l.p1.x << "," << l.p1.y << ") -> ("
                  << l.p2.x << "," << l.p2.y << ")\n";
    }

    if (options.svgFile) {
        exportToSvg(contour, lines, *options.svgFile);
        std::cerr << "SVG written to: " << *options.svgFile << "\n";
    }

    return 0;
}
