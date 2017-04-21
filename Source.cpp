#include <cv.h>
#include <highgui.h>
#include <iostream>
#include<opencv2/opencv.hpp>
#include<stdio.h>
#include <stdlib.h>  
#include <time.h>
#include<math.h>
#include<fstream>


using namespace std;
using namespace cv;

double gRemoveSize, gBigObj;
int gProtectSize = 51;
string imgName = "08430000100A1083";

const double RANGE_BEGIN = 0.0;
const double RANGE_END = 180.0;
const double MINUS = 0.1;

#define PI 3.1416
#define THRESHOLDofSIZE 5000 //面積多少以上要保留
#define ThetaNum 180 //0~180窮舉
#define THRESHOLD 50
#define MAX_OF_COUNT 2 //角度數量前幾大的角度要被計算
#define THRESHOLD_SETTING 20
//#define RANGE_BEGIN 0.0 //角度下限
//#define RANGE_END 180.0 //角度上限
//#define MINUS 0.1 //角度間隔

void calcLinesP(const Mat &input, std::vector<Vec4i> &lines);
void drawLinesP(Mat &input, const std::vector<Vec4i> &lines);
void calcLines(const Mat &input, std::vector<Vec2f> &lines);
void drawLines(Mat &input, const std::vector<Vec2f> &lines);

//切割圖片1/10
Mat CutInput(Mat input){
	Mat dst = input(Rect(input.cols - input.cols / 20, 0, input.cols / 20 - 1, input.rows));
	return dst;
}

int main(){
	string temp = imgName + ".png";

	Mat img = imread(temp, 0);
	threshold(img, img, 128, 255, THRESH_BINARY);
	img = CutInput(img);
	img = ~img;
	imwrite("(0)切割原圖.png", img);

	Mat white(img.rows, img.cols, CV_8UC1, Scalar(255));

	//去除小物件
	Mat labelImage, stats, centroids;
	int nLabels = connectedComponentsWithStats(img, labelImage, stats, centroids, 8, CV_32S);
	int label;
	vector<int>labelForEliminate;

	for (int label = 1; label < nLabels; ++label){
		if (stats.at<int>(label, cv::CC_STAT_AREA) > THRESHOLDofSIZE){
			labelForEliminate.push_back(label);
		}
	}

	for (int i = 0; i < img.rows; i++)
	for (int j = 0; j < img.cols; j++){
		for (int k = 0; k < labelForEliminate.size(); k++){
			if (labelImage.at<int>(i, j) == labelForEliminate[k]){
				white.at<uchar>(i, j) = 0;
			}
		}
	}

	imwrite("(1)去除小物件圖.png", white);

	//細化--
	Mat tempdst;

	int Z_count, T_count, flag = 100;

	threshold(white, white, 125, 255, THRESH_BINARY);//二值化
	copyMakeBorder(white, tempdst, 1, 1, 1, 1, BORDER_REPLICATE);//鏡像
	//imshow("input", white);
	//waitKey();

	while (flag != 0){
		for (int i = 1; i < white.rows + 1; i++)
		for (int j = 1; j < white.cols + 1; j++){
			if ((int)tempdst.at<uchar>(i, j) == 0){
				Z_count = 0;
				T_count = 0;
				int neighbor[8];
				neighbor[0] = (int)tempdst.at<uchar>(i - 1, j - 1);
				neighbor[1] = (int)tempdst.at<uchar>(i - 1, j);
				neighbor[2] = (int)tempdst.at<uchar>(i - 1, j + 1);
				neighbor[3] = (int)tempdst.at<uchar>(i, j + 1);
				neighbor[4] = (int)tempdst.at<uchar>(i + 1, j + 1);
				neighbor[5] = (int)tempdst.at<uchar>(i + 1, j);
				neighbor[6] = (int)tempdst.at<uchar>(i + 1, j - 1);
				neighbor[7] = (int)tempdst.at<uchar>(i, j - 1);

				for (int k = 0; k < 8; k++){
					if (neighbor[k] < 125) Z_count++;
					if (k != 0){//1~8
						if (neighbor[k] - neighbor[k - 1] >= 125) T_count++;
					}
					if (k == 0)if (neighbor[0] - neighbor[7] >= 125)T_count++;
				}
				if (Z_count <= 6 && Z_count >= 2 && T_count == 1){
					if (flag % 2 == 1 && (neighbor[1] == 255 || neighbor[3] == 255 || neighbor[5] == 255) && (neighbor[1] == 255 || neighbor[7] == 255 || neighbor[3] == 255)){
						white.at<uchar>(i - 1, j - 1) = 255;
						//tempdst.at<uchar>(i, j) = 255;

					}
					if (flag % 2 == 0 && (neighbor[7] == 255 || neighbor[3] == 255 || neighbor[5] == 255) && (neighbor[1] == 255 || neighbor[7] == 255 || neighbor[5] == 255)){
						white.at<uchar>(i - 1, j - 1) = 255;
						//tempdst.at<uchar>(i, j) = 255;
					}
				}
			}
		}//for
		/*system("cls");
		cout << 101 - flag << "%";*/

		copyMakeBorder(white, tempdst, 1, 1, 1, 1, BORDER_REPLICATE);
		flag--;
	}//while
	imwrite("(2)細化.png", white);

	//---


	//手寫測線
	Mat dst=white.clone();
	//Canny(white, dst, 50, 150);
	//imwrite("canny.png", dst);
	//dst = ~dst;

	//法距長範圍0~整張圖最長邊*2(值角三角最長邊) ?這裡用int可以嗎待確認? *****(乘以二因為套霍夫公式會算出負值 要把她右移變正)*****
	int DistNum = 2 * sqrt(dst.rows*dst.rows + dst.cols*dst.cols);
	//投票箱製作 y軸ThetaNum x軸為法距長 
	int **box;
	box = new int*[(int)((RANGE_END - RANGE_BEGIN) / MINUS)];

	for (double i = 0; i < RANGE_END - RANGE_BEGIN; i += MINUS){
		box[(int)(i*(1 / MINUS))] = new int[DistNum];
	}
	//投票箱歸零
	for (double i = 0; i < RANGE_END - RANGE_BEGIN; i += MINUS)
	for (int j = 0; j < DistNum; j++){
		box[(int)(i*(1 / MINUS))][j] = 0;
	}
	//轉換極作標
	//整張圖遍歷
	for (int i = 0; i < dst.rows; i++)
	for (int j = 0; j < dst.cols; j++){
		//只對黑點有興趣
		if ((int)dst.at<uchar>(i, j) == 0){
			//0~179代入
			for (double x = RANGE_BEGIN; x < RANGE_END; x += MINUS){
				int r;
				//核心公式
				r = round(i*cos((double)x*(PI / 180)) + (double)j*sin(x*(PI / 180))) + DistNum / 2;//r錯了 有負值 所以要加DistNum/2右移
				box[(int)((x - RANGE_BEGIN)*(1.0 / MINUS))][r] += 1;
			}
		}
	}
	//每個角度數量累計
	ofstream outfile_2("箱子角度統計.txt");
	double MAXangle = -1, tempMax = -1;
	int *angle_array;
	angle_array = new int[(int)((RANGE_END - RANGE_BEGIN) / MINUS)];
	for (int i = 0; i < (int)((RANGE_END - RANGE_BEGIN) / MINUS); i++)
		angle_array[i] = 0;
	//int angle_array[(int)((RANGE_END - RANGE_BEGIN) / MINUS)] = { 0 };//0~180度每個角度的數量

	for (double i = 0; i < RANGE_END - RANGE_BEGIN; i += MINUS){
		for (int j = 0; j < DistNum; j++){
			if (box[(int)(i*(1.0 / MINUS))][j] > THRESHOLD){
				angle_array[(int)(i*(1.0 / MINUS))] = angle_array[(int)(i*(1.0 / MINUS))] + box[(int)(i*(1.0 / MINUS))][j];
			}
		}
		cout << (int)(i*(1 / MINUS)) << "=" << (double)(i + RANGE_BEGIN) << endl;
		if (tempMax < angle_array[(int)(i*(1.0 / MINUS))] && i*(1.0 / MINUS) != (90 - RANGE_BEGIN)*(1.0 / MINUS)){//不記九十度 因為90太多了
			MAXangle = i;
			tempMax = angle_array[(int)(i*(1.0 / MINUS))];
		}
		//outfile_2 << i + RANGE_BEGIN << "度=" << angle_array[(int)(i*(1 / MINUS))] << "個" << endl;
	}
	//cout << "max=" << MAXangle + RANGE_BEGIN << endl;


	//畫(new)----
	double maxcase = -1, maxANGLE = -1, maxR = -1;
	Mat testoutput(img.rows, img.cols, CV_8UC3, Scalar(255, 255, 255));
	for (int i = 0; i < img.rows; i++)
	for (int j = 0; j < img.cols; j++){
		//if ((int)dst.at<uchar>(i, j) == 0)testoutput.at < Vec3b >(i, j) = 0;
	}
	//找最大票箱
	for (double i = 0; i < RANGE_END - RANGE_BEGIN; i += MINUS){
		for (int j = 0; j < DistNum; j++){
			if (box[(int)(i*(1.0 / MINUS))][j]>maxcase){
				maxcase = box[(int)(i*(1.0 / MINUS))][j];
				maxANGLE = i;
				maxR = j;
			}
		}
	}
	cout << "最大度數=" << maxANGLE + RANGE_BEGIN << endl;

	for (int i = 0; i < img.rows; i++)
	for (int j = 0; j < img.cols; j++){
		int r;
		//核心公式
		r = round(i*cos((double)maxANGLE*(PI / 180)) + (double)j*sin(maxANGLE*(PI / 180))) + DistNum / 2;//r錯了 有負值 所以要加DistNum/2右移
		if (r == maxR)testoutput.at<Vec3b>(i, j) = 0;
	}
	imwrite("(3)畫出測到最長線.png", testoutput);

	

	//轉圖
	Mat src = imread(temp, 1);
	Mat dst2(src.rows, src.cols,CV_8UC1,Scalar(255));

	//選定幾何轉換前後相對的三個點
	Point2f srcTri[3];
	srcTri[0] = Point2f(0, 0);
	srcTri[1] = Point2f(src.cols - 1, 0);
	srcTri[2] = Point2f(0, src.rows - 1);

	Point2f dstTri[3];
	dstTri[0] = Point2f(0, src.rows*0.3);
	dstTri[1] = Point2f(src.cols*0.8, 0);
	dstTri[2] = Point2f(src.cols*0.1, src.rows*0.9);

	cout << (maxANGLE + RANGE_BEGIN - 90);
	//設定旋轉中心、旋轉角度和縮放倍率
	Point center = Point(dst2.cols / 2, dst2.rows / 2);
	double angle = -(maxANGLE + RANGE_BEGIN - 90);
	double scale = 1.0;

	Mat rot_mat = getRotationMatrix2D(center, angle, scale);
	warpAffine(src, dst2, rot_mat, dst2.size(), INTER_LINEAR, BORDER_CONSTANT, Scalar(255,255,255));

	//imshow("origin", src);
	//imshow("Affine_1", dst1);
	//imshow("Affine_2", dst2);
	imwrite("(4)轉完圖結果.png", dst2);

	/*//opencv函式呼叫
	Mat result1 = white.clone();
	Mat result2 = white.clone();

	vector<Vec2f> lines;
	calcLines(white, lines);
	drawLines(result2, lines);
	vector<Vec4i> linesP;
	calcLinesP(white, linesP);
	drawLinesP(result1, linesP);

	imwrite("機率霍夫結果.png", result1);
	imwrite("一般霍夫結果.png", result2);
	*/

	system("pause");
	return 0;
}

void calcLinesP(const Mat &input, std::vector<Vec4i> &lines){
	Mat contours;
	Canny(input, contours, 50, 150);
	imwrite("canny.png", contours);
	lines.clear();
	HoughLinesP(contours, lines,1, CV_PI / 180, input.rows / 4);
}


void drawLinesP(Mat &input, const std::vector<Vec4i> &lines){
	for (int i = 0; i<lines.size(); i++){
		line(input, Point(lines[i][0], lines[i][3]), Point(lines[i][4], lines[i][5]), Scalar(255, 0, 0), 3);
	}
}

void calcLines(const Mat &input, std::vector<Vec2f> &lines){
	Mat contours;
	Canny(input, contours, 50, 150);
	lines.clear();
	HoughLines(contours, lines, 1, CV_PI / 180, input.rows / 4);
}

void drawLines(Mat &input, const std::vector<Vec2f> &lines){
	for (int i = 0; i<lines.size(); i++){
		double r = lines[i][0];
		double theta = lines[i][6];
		if (theta<PI / 4.0 || theta>3 * PI / 4.0){
			Point pt1(r / cos(theta), 0);
			Point pt2((r - input.rows*sin(theta)) / cos(theta), input.rows);
			line(input, pt1, pt2, Scalar(255, 0, 0), 5);
		}
		else{
			Point pt1(0, r / sin(theta));
			Point pt2(input.cols, (r - input.cols*cos(theta)) / sin(theta));
			line(input, pt1, pt2, Scalar(255, 0, 0), 3);
		}
	}
}
