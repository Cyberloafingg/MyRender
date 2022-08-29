#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include <iostream>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model* model = NULL;
const int width = 800;
const int height = 800;

void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
    bool steep = false;
    if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1 - x0;
    int dy = y1 - y0;
    int derror2 = std::abs(dy) * 2;
    int error2 = 0;
    int y = y0;
    for (int x = x0; x <= x1; x++) {
        if (steep) {
            image.set(y, x, color);
        }
        else {
            image.set(x, y, color);
        }
        error2 += derror2;
        if (error2 > dx) {
            y += (y1 > y0 ? 1 : -1);
            error2 -= dx * 2;
        }
    }
}

//void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
//    // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!) 
//    if (t0.y > t1.y) std::swap(t0, t1);
//    if (t0.y > t2.y) std::swap(t0, t2);
//    if (t1.y > t2.y) std::swap(t1, t2);
//    int total_height = t2.y - t0.y;
//    for (int y = t0.y; y <= t1.y; y++) {
//        int segment_height = t1.y - t0.y + 1;
//        float alpha = (float)(y - t0.y) / total_height;
//        float beta = (float)(y - t0.y) / segment_height; // be careful with divisions by zero 
//        Vec2i A = t0 + (t2 - t0) * alpha;
//        Vec2i B = t0 + (t1 - t0) * beta;
//        if (A.x > B.x) std::swap(A, B);
//        for (int j = A.x; j <= B.x; j++) {
//            image.set(j, y, color); // attention, due to int casts t0.y+i != A.y 
//        }
//    }
//    for (int y = t1.y; y <= t2.y; y++) {
//        int segment_height = t2.y - t1.y + 1;
//        float alpha = (float)(y - t0.y) / total_height;
//        float beta = (float)(y - t1.y) / segment_height; // be careful with divisions by zero 
//        Vec2i A = t0 + (t2 - t0) * alpha;
//        Vec2i B = t1 + (t2 - t1) * beta;
//        if (A.x > B.x) std::swap(A, B);
//        for (int j = A.x; j <= B.x; j++) {
//            image.set(j, y, color); // attention, due to int casts t0.y+i != A.y 
//        }
//    }
//    //for (int i = 0; i < total_height; i++) {
//    //    bool half = i > t1.y - t0.y || t1.y == t0.y;
//    //    int segment_height = half ? t2.y - t1.y : t1.y - t0.y;
//    //    float alpha = (float)i / total_height;
//    //    float beta = (float)(i - (half ? t1.y - t0.y : 0)) / segment_height;
//    //    Vec2i A = t0 + (t2 - t0) * alpha;
//    //    Vec2i B = half ? t1 + (t2 - t1) * beta : t0 + (t1 - t0) * beta;
//    //    if (A.x > B.x) std::swap(A, B);
//    //    for (int j = A.x; j <= B.x; j++) {
//    //        image.set(j, t0.y + i, color); // attention, due to int casts t0.y+i != A.y 
//    //    }
//    //}
//
//}



void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
    Vec2i t[3] = { t0, t1, t2 };
    Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
    Vec2i bboxmax(0, 0);
    for (int i = 0; i < 3; i++) {
        bboxmin.x = std::max(0, std::min(bboxmin.x, t[i].x));
        bboxmin.y = std::max(0, std::min(bboxmin.y, t[i].y));

        bboxmax.x = std::min(image.get_width() - 1, std::max(bboxmax.x, t[i].x));
        bboxmax.y = std::min(image.get_height() - 1, std::max(bboxmax.y, t[i].y));
    }
    // 遍历box区域，判断像素是否在三角内
    Vec2i P,AP,BP,CP;
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
            AP = P - t0; BP = P - t1; CP = P - t2;
            Vec2i arr[3] = {AP, BP, CP};
            if( ((arr[0].x * arr[1].y - arr[0].y * arr[1].x)*(arr[1].x * arr[2].y - arr[1].y * arr[2].x)) >= 0 && 
                ((arr[2].x * arr[0].y - arr[2].y * arr[0].x)*(arr[0].x * arr[1].y - arr[0].y * arr[1].x)) >= 0){
                image.set(P.x, P.y, color);
            }
            //int ts = 1;
            //for (int i = 0; i < 3; i++) {
            //    ts *= (arr[i].x * arr[(i + i) % 3].y - arr[i].y * arr[(i + 1) % 3].x);
            //}
            //if (ts > 0) {
            //    image.set(P.x, P.y, color);
            //}
        }
    }
}




int main(int argc, char** argv) {
    if (2 == argc) {
        model = new Model(argv[1]);
    }
    else {
        model = new Model("F:/CGraphic/TinyRender/TinyRender/Project1/obj/african_head.obj");
    }
    // 绘制模型
    TGAImage image(width, height, TGAImage::RGB);

    Vec3f light_dir(0, 0, -1);
    bool flag;
    int count = 0;
    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3];
        Vec3f world_coords[3];
        for (int j = 0; j < 3; j++) {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
            world_coords[j] = v;
        }
        // 求叉乘
        std::cout << i << std::endl;
        std::cout << "A:" << world_coords[0] << "B:" << world_coords[1] << "C:" << world_coords[2];
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        std::cout << "Z:"  << n.z << std::endl;
        if (n.z > 0) count++;
        n.normalize();
        float intensity = n * light_dir;
        /*std::cout << intensity << std::endl << std::endl;*/
        if (n.z > 0) {
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
            if (i == (model->nfaces()-2) || i == 993 || i == 995) {
                triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, red);
            }
        }
    }
    std::cout << count << std::endl;
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}
    