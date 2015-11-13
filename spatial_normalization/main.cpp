//�̑���	��ԓI�W����
//�쐬��	������y�H�C���c��y�H
//�@����	�~��
/*�擪�œ��{���ł�����ł����΃\�[�X�c���[�ŕ\�������Ƃ��ɕ����������Ȃ��炵���̂�*/
#include "naricommon.h"
#include "ip/narigaussian.h"

#include "naritimer.h"
#include "naricaseinfo.h"
#include <algorithm>
#include "other/narimhd.h"

#include <sstream>
#include "naricommon.h"
#include "narisystem.h"
#include "naripathline.h"
#include "naricaseinfo.h"
#include "ip/narimorphology.h"
#include "ip/narilabeling.h"

#include "ip/naricontour.h"
#include "ip/naridistance.h"
#include "other/narimhd.h"
#include "naritimer.h"
#include "ip/narirbf.h"
#include "../mist/vector.h"
#include "narivectorpp.h"
#include "ip/nariinterpolate.h"

#include "raw_io.h"
#include "info.h"
#include "template_matching.h"


template <class IMG_T>
void deformation_using_movement(int xeRef, int yeRef, int zeRef, int xeFl, int yeFl, int zeFl, nari::vector<float> &imgMoveX, nari::vector<float> &imgMoveY, nari::vector<float> &imgMoveZ, nari::vector<IMG_T> &imgI, nari::vector<IMG_T> &imgO)
{
	for(int z = 0; z < zeRef; z++)
	{
		for(int y = 0; y < yeRef; y++)
		{
			for(int x = 0; x < xeRef; x++)
			{
				int s = z * xeRef * yeRef + y * xeRef + x;
				imgO[s] = static_cast<IMG_T>(nari::interpolate_value::linear(imgI.ptr(), imgMoveX[s], imgMoveY[s], imgMoveZ[s], xeFl, yeFl, zeFl));
			}
		}
	}
}
template <class LABEL_T>
void voxelization(nari::mhd &mhdLabel, nari::vector<LABEL_T> &imgLabel)
{

	// voxelization
	int xe = mhdLabel.size1();
	int ye = mhdLabel.size2();
	int ze = mhdLabel.size3();

	int xe_o = xe;
	int ye_o = ye;
	int ze_o = ze;
	double xr = mhdLabel.reso1();
	double yr = mhdLabel.reso2();
	double zr = mhdLabel.reso3();

	double zr_r = xr;
	int ze_r = (int)(ze * zr / zr_r);

	// voxelize
	std::cout << "voxelization : " << ze << " > " << ze_r << ", " << zr << " > " << zr_r << std::endl;
	ze = ze_r;
	zr = zr_r;

	// shrink
	int sh = 2;
	std::cout << "shrink : 1/" << sh << std::endl;
	int xe_r = xe / sh;
	int ye_r = ye / sh;
	ze_r = ze / sh;

	double xr_r = xr * sh;
	double yr_r = yr * sh;
	zr_r = zr * sh;
	nari::vector< float > imgDist(xe_o * ye_o * ze_o);
	nari::distance::euclidean_signed_distance_transform(imgLabel.ptr(), imgDist.ptr(), xe_o, ye_o, ze_o);
	nari::vector< float > img_shrink;
	nari::interpolate::linear(imgDist, img_shrink, xe_o, ye_o, ze_o, xe_r, ye_r, ze_r);
	imgLabel.resize(xe_r * ye_r * ze_r);
	for(int s = 0; s < imgLabel.size(); s++)
	{
		if(img_shrink[s] < 0) imgLabel[s] = 1;
		else imgLabel[s] = 0;
	}
	mhdLabel.size123(xe_r, ye_r, ze_r);
	mhdLabel.reso123(xr_r, yr_r, zr_r);
}

void main(int argc, char *argv[])
{
	//�����炭rbf�ό`�Ɏg���^�̒�`�C3�����p���H
	typedef nari::rbf<double, double, 3> rbf_t;
	typedef rbf_t::disp_t disp_t_3;
	typedef rbf_t::disp_list_t disp_list_t_3;


	//�e�L�X�g��path���e���v���[�g�}�b�`���O�̃p�����[�^���ǂ߂�悤�ɂ��܂���
	if (argc != 2)
	{
		std::cout << "usage : make_filter.exe input_info.txt" << std::endl << "Hit enter-key to exit...";
		getchar();
		exit(1);
	} 

	text_info input_info;
	input_info.input(argv[1]);



	// information of case
	nari::case_list_t cases;
	cases.load_file_txt(input_info.pathId, true);


	//��Ǘ�̏����擾(�T�C�Y�C����_�Ȃ�)
	//nari::case_t caseRef("jamit09", "jamit09_08");
	nari::mhd mhdRef;
	disp_list_t_3 listDisp;
	mhdRef.load( input_info.dirBase + input_info.caseRef.dir + "/dispImg/" + input_info.caseRef.basename + ".mhd");


	int xeRef = mhdRef.size1();
	int yeRef = mhdRef.size2();
	int zeRef = mhdRef.size3();
	double xrRef = mhdRef.reso1();
	double yrRef = mhdRef.reso2();
	double zrRef = mhdRef.reso3();
	const double shrinkValue = 3.0;	//1/3�ɏk�����܂�
	int xeRefS = (int)(xeRef / shrinkValue + 0.5);
	int yeRefS = (int)(yeRef / shrinkValue + 0.5);
	int zeRefS = (int)(zeRef / shrinkValue + 0.5);
	double xrRefS = xrRef * (double)xeRef / xeRefS;
	double yrRefS = yrRef * (double)yeRef / yeRefS;
	double zrRefS = zrRef * (double)zeRef / zeRefS;

	nari::mhd mhdRefS;
	mhdRefS = mhdRef;
	mhdRefS.size123(xeRefS, yeRefS, zeRefS);
	mhdRefS.reso123(xrRefS, yrRefS, zrRefS);

	//���������Ǘ�̃e���v���[�g�}�b�`���O���s���v���_�̍��W���Z�o,�@�ۑ�
	std::ofstream Ref_list(input_info.dirBase +input_info. caseRef.dir + "/disp/" + input_info.caseRef.basename + ".txt");

	nari::vector<int> x_disp (xeRef,0);
	nari::vector<int> y_disp (yeRef,0);
	nari::vector<int> z_disp (zeRef,0);

	int p,q,r;
	int x,y,z;
	p=0;
	q=0;
	r=0;
	int rx = input_info.rangex;
	int ry = input_info.rangey;
	int rz = input_info.rangez;

	for ( x = rx ; x < xeRef ; x+=rx ){
		x_disp[p] = x;
		p++;
	} 

	for ( y = ry ; y < yeRef ; y+=ry ){
		y_disp[q] = y;
		q++;
	}  

	for ( z = rz ; z < zeRef ; z+=rz ){
		z_disp[r] = z;
		r++;
	}   

	for ( x = 0; x <= p; x++ ){
		for ( y = 0; y <= q; y++ ){
			for ( z = 0; z <= r; z++ ){
				Ref_list << x_disp[x]  << std::endl;
				Ref_list << y_disp[y]  << std::endl;
				Ref_list << z_disp[z]  << std::endl;
			}
		}
	}

	//�v���_�̍��W�����[�h,listDisp�ɑ��
	rbf_t::load_cordinates(input_info.dirBase +input_info. caseRef.dir + "/disp/" + input_info.caseRef.basename + ".txt", listDisp);

	
	//�Ǘ჋�[�v
	for(int i = 0; i < cases.size(); i++)
	{
		std::cout << "case : " << i + 1 << std::endl;


		//�����Ǘ�̏����擾(�T�C�Y�C����_�Ȃ�)
		nari::mhd mhdFl;
		mhdFl.load( input_info.dirOrg + cases[i].dir + cases[i].basename + "_3.mhd");
		int xeFl = mhdFl.size1();
		int yeFl = mhdFl.size2();
		int zeFl = mhdFl.size3();
		double xrFl = mhdFl.reso1();
		double yrFl = mhdFl.reso2();
		double zrFl = mhdFl.reso3();
		int xeFlS = (int)(xeFl / shrinkValue + 0.5);
		int yeFlS = (int)(yeFl / shrinkValue + 0.5);
		int zeFlS = (int)(zeFl / shrinkValue + 0.5);
		double xrFlS = xrFl * (double)xeFl / xeFlS;
		double yrFlS = yrFl * (double)yeFl / yeFlS;
		double zrFlS = zrFl * (double)zeFl / zeFlS;

		nari::mhd mhdFlS;
		mhdFlS = mhdFl;
		mhdFlS.size123(xeFlS, yeFlS, zeFlS);
		mhdFlS.reso123(xrFlS, yrFlS, zrFlS);


		//��Ǘ�̓ǂݍ���
		nari::vector<short> imgRef(xeRef * yeRef * zeRef);
		mhdRef.save_mhd_and_image(imgRef, input_info.dirO + cases[i].dir + "org/" + cases[i].basename + "_3" + ".raw");
		//�����Ǘ�̓ǂݍ���
		nari::vector<short> imgFl(xeFl * yeFl * zeFl);
		mhdFl.load_mhd_and_image(imgFl, input_info.dirOrg + cases[i].dir + cases[i].basename + "_3" + ".raw");
		//�����Ǘ�̌v���_���i�[����x�N�g����p��
		nari::vector<int> x_disp2 (xeRef,0);
		nari::vector<int> y_disp2 (yeRef,0);
		nari::vector<int> z_disp2 (zeRef,0);

		//�e���v���[�g�}�b�`���O
		template_mathcing(imgRef, imgFl ,x_disp ,y_disp ,z_disp ,x_disp2 ,y_disp2 ,z_disp2 ,xeRef,yeRef, zeRef,
						 xeFl, yeFl, zeFl, input_info.tmp , rx ,ry, rz);
		//�e���v���[�g�}�b�`���O�ɂ�茈�肵���Ή��_�̍��W���e�L�X�g�ɕۑ�
		std::ofstream Fl_list(input_info.dirDisp + cases[i].dir + "disp/" + cases[i].basename + ".txt");

		for ( x = 0; x <= p; x++ ){
			for ( y = 0; y <= q; y++ ){
				for ( z = 0; z <= r; z++ ){
					Fl_list << x_disp2[x]  << std::endl;
					Fl_list << y_disp2[y]  << std::endl;
					Fl_list << z_disp2[z]  << std::endl;
				}
			}
		}
		//�e���v���[�g�}�b�`���O�̕����Ǘ�̌v���_�����[�h
		rbf_t::load_values(input_info.dirDisp + cases[i].dir + "disp/" + cases[i].basename + ".txt", listDisp);	//����_�̃��[�h

		mist::matrix<double> wMat = rbf_t::get_w_mat(listDisp, nari::rbf_kernel::biharmonic());

		//�ړ���i�[�p�i��Ǘ�̍��W����Ή����镂���Ǘ�̍��W���擾�j
		nari::vector<float> imgMoveX(xeRefS * yeRefS * zeRefS);
		nari::vector<float> imgMoveY(xeRefS * yeRefS * zeRefS);
		nari::vector<float> imgMoveZ(xeRefS * yeRefS * zeRefS);

		{
			int xe = xeRefS;
			int ye = yeRefS;
			int ze = zeRefS;
			double sX = (double)xeRef / xe;
			double sY = (double)yeRef / ye;
			double sZ = (double)zeRef / ze;
			std::cout << "num Displacement : " << listDisp.size() << std::endl;
			for(int z = 0; z <ze; z++)
			{
				printf("z : [%5d/%5d]\r", z, ze - 1);
				for(int y = 0; y < ye; y++)
				{
					for(int x = 0; x < xe; x++)
					{
						int s = z * xe * ye + y * xe + x;
						mist::matrix<double> cordinate = rbf_t::rbf_interpolation(listDisp, wMat, x * sX, y * sY, z * sZ, nari::rbf_kernel::biharmonic());
						imgMoveX[s] = (float)cordinate[0];
						imgMoveY[s] = (float)cordinate[1];
						imgMoveZ[s] = (float)cordinate[2];
					}
				}
			}
		}
		std::cout << std::endl;
		//���ȉ��ɓ_���Ƃ̈ړ��������W���i�[����
		mhdRefS.save_mhd_and_image(imgMoveX, input_info.dirO + cases[i].dir + "move/" + "x/" + cases[i].basename + ".raw");
		mhdRefS.save_mhd_and_image(imgMoveY, input_info.dirO + cases[i].dir + "move/" + "y/" + cases[i].basename + ".raw");
		mhdRefS.save_mhd_and_image(imgMoveZ, input_info.dirO + cases[i].dir + "move/" + "z/" + cases[i].basename + ".raw");

		//���`�ϊ��Cnari::interpolate::linear�ł�nari�^��vector����mist�^��array�������ɕϊ�,mist�̐��`��Ԃ��g�p���Ă���
		//xeRef�͌��̃T�C�Y�Cxerefs��1/3�̑傫���C1/3�̑傫���̉摜�����̃T�C�Y�ɖ߂��Ă���
		nari::interpolate::linear(imgMoveX, imgMoveX, xeRefS, yeRefS, zeRefS, xeRef, yeRef, zeRef);
		nari::interpolate::linear(imgMoveY, imgMoveY, xeRefS, yeRefS, zeRefS, xeRef, yeRef, zeRef);
		nari::interpolate::linear(imgMoveZ, imgMoveZ, xeRefS, yeRefS, zeRefS, xeRef, yeRef, zeRef);


		nari::vector<short> imgI(xeFl * yeFl * zeFl), imgO(xeRef * yeRef * zeRef);
		mhdFl.load_mhd_and_image(imgI, input_info.dirOrg + cases[i].dir + cases[i].basename + "_3" + ".mhd");

		deformation_using_movement(xeRef, yeRef, zeRef, xeFl, yeFl, zeFl, imgMoveX, imgMoveY, imgMoveZ, imgI, imgO);
		mhdRef.save_mhd_and_image(imgO, input_info.dirO + cases[i].dir + "org/" + cases[i].basename + "_3" + ".raw");

		const	int nLabel = 8;
		const	std::string	nameLabel[] = {"liver/", "heart/", "kidney/", "gallbladder/", "stomach/", "tumor/", "l_kidney/", "spleen/"};
		//���탉�x�������Ԃɕό`���Ă���
		for(int l = 0; l < nLabel; l++)
		{
			nari::vector<unsigned char> imgLabel;
			nari::mhd mhdLabel;
			std::ostringstream oss;
			oss << input_info.dirLabel << cases[i].dir + nameLabel[l] + cases[i].basename << "_3.mhd";

			if( !nari::system::file_is_exist(oss.str()) ) continue;
			mhdLabel.load_mhd_and_image(imgLabel, oss.str());
			voxelization(mhdLabel, imgLabel);
			mhdLabel.size123(xeFl, yeFl, zeFl);
			mhdLabel.save_mhd_and_image(imgLabel, input_info.dirO + cases[i].dir + "cut/" + nameLabel[l] + cases[i].basename + ".raw");

			// �����ϊ�
			for( int s = 0; s < xeFl * yeFl; s++)	imgLabel[s] = 0;
			for( int s = xeFl * yeFl * ( zeFl - 1); s < xeFl * yeFl * zeFl; s++)	imgLabel[s] = 0;
			nari::vector<float> imgLabelDist(xeFl * yeFl * zeFl), imgLabelDistO(xeRef * yeRef * zeRef);
			nari::distance::euclidean_signed_distance_transform(imgLabel.ptr(), imgLabelDist.ptr(), xeFl, yeFl, zeFl);

			//�ό`
			deformation_using_movement(xeRef, yeRef, zeRef, xeFl, yeFl, zeFl, imgMoveX, imgMoveY, imgMoveZ, imgLabelDist, imgLabelDistO);

			oss.str("");

			//�����摜��臒l����
			nari::vector<unsigned char> imgLabelO(xeRef * yeRef * zeRef, 0);
			for(int s = 0; s < xeRef * yeRef * zeRef; s++)
			{
				if(imgLabelDistO[s] < 0) imgLabelO[s] = 1;
			}
			oss.str("");

			//���x���̕ۑ�
			oss << input_info.dirO << cases[i].dir + nameLabel[l] + cases[i].basename << ".raw";
			mhdRef.save_mhd_and_image(imgLabelO, oss.str());

		}
	}
}