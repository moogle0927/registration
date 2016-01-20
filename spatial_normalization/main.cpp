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
	for (int z = 0; z < zeRef; z++)
	{
		for (int y = 0; y < yeRef; y++)
		{
			for (int x = 0; x < xeRef; x++)
			{
				int s = z * xeRef * yeRef + y * xeRef + x;
				imgO[s] = static_cast<IMG_T>(nari::interpolate_value::linear(imgI.ptr(), imgMoveX[s], imgMoveY[s], imgMoveZ[s], xeFl, yeFl, zeFl));
			}
		}
	}
}
template <class LABEL_T>

//���̏����͂Ƃ肠�����v��Ȃ�
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
	//�������������s�������ze��ze_r�Ƃ���
	int ze_r = (int)(ze * zr / zr_r);

	// voxelize
	std::cout << "voxelization : " << ze << " > " << ze_r << ", " << zr << " > " << zr_r << std::endl;
	ze = ze_r;
	zr = zr_r;

	// shrink
	//�摜�̑傫����1/sh�ɏk��
	int sh = 2;
	std::cout << "shrink : 1/" << sh << std::endl;
	int xe_r = xe / sh;
	int ye_r = ye / sh;
	ze_r = ze / sh;

	double xr_r = xr * sh;
	double yr_r = yr * sh;
	zr_r = zr * sh;
	nari::vector< float > imgDist(xe_o * ye_o * ze_o);
	//���͂������x���摜���狗���ϊ��摜���쐬
	nari::distance::euclidean_signed_distance_transform(imgLabel.ptr(), imgDist.ptr(), xe_o, ye_o, ze_o);
	nari::vector< float > img_shrink;
	//��������1/sh�ɂ��������ϊ��摜���쐬
	nari::interpolate::linear(imgDist, img_shrink, xe_o, ye_o, ze_o, xe_r, ye_r, ze_r);
	imgLabel.resize(xe_r * ye_r * ze_r);
	//�����ϊ��摜��臒l�������ă��x���摜�ɖ߂�
	for (int s = 0; s < imgLabel.size(); s++)
	{
		if (img_shrink[s] < 0) imgLabel[s] = 1;
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
	nari::mhd mhdRef;
	disp_list_t_3 listDisp;
	mhdRef.load(input_info.dirRef + input_info.caseRef_dir + input_info.caseRef_name + ".mhd");

	int xeRef = mhdRef.size1();
	int yeRef = mhdRef.size2();
	int zeRef = mhdRef.size3();
	double xrRef = mhdRef.reso1();
	double yrRef = mhdRef.reso2();
	double zrRef = mhdRef.reso3();
		
	//�v���_�̍��W�����[�h,listDisp�ɑ��
	rbf_t::load_cordinates(input_info.dirRef + input_info.caseRef_dir + "disp/" + input_info.caseRef_name + ".txt", listDisp);

	//�Ǘ჋�[�v
	for (int i = 0; i < cases.size(); i++)
	{
		std::cout << "case : " << i + 1 << std::endl;

		//�����Ǘ�̏����擾(�T�C�Y�C����_�Ȃ�)
		nari::mhd mhdFl;
		mhdFl.load(input_info.dirFl + cases[i].dir + cases[i].basename + ".mhd");


		int xeFl = mhdFl.size1();
		int yeFl = mhdFl.size2();
		int zeFl = mhdFl.size3();
		double xrFl = mhdFl.reso1();
		double yrFl = mhdFl.reso2();
		double zrFl = mhdFl.reso3();

		std::cout << "(^^)<�摜�ǂݍ��ނ�" << std::endl;

		//�����Ǘ�̌v���_�����[�h
		rbf_t::load_values(input_info.dirFl + cases[i].dir + "disp/" + cases[i].basename + ".txt", listDisp);	//����_�̃��[�h
		std::cout << "���[�h����" << std::endl;
		mist::matrix<double> wMat = rbf_t::get_w_mat(listDisp, nari::rbf_kernel::biharmonic(), input_info.lamb);
		std::cout << wMat << std::endl;

		//�ړ���i�[�p�i��Ǘ�̍��W����Ή����镂���Ǘ�̍��W���擾�j
		nari::vector<float> imgMoveX(xeRef * yeRef * zeRef);
		nari::vector<float> imgMoveY(xeRef * yeRef * zeRef);
		nari::vector<float> imgMoveZ(xeRef * yeRef * zeRef);

		{
			int xe = xeRef;
			int ye = yeRef;
			int ze = zeRef;
			double sX = (double)xeRef / xe;
			double sY = (double)yeRef / ye;
			double sZ = (double)zeRef / ze;
			std::cout << "num Displacement : " << listDisp.size() << std::endl;
			for (int z = 0; z < ze; z++)
			{
				printf("z : [%5d/%5d]\r", z, ze - 1);
				for (int y = 0; y < ye; y++)
				{
					for (int x = 0; x < xe; x++)
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
		
		//���ȉ��ɓ_���Ƃ̈ړ��������W���i�[����
		mhdRef.save_mhd_and_image(imgMoveX, input_info.dirO + cases[i].dir + "move/" + "x/" + cases[i].basename + ".raw");
		mhdRef.save_mhd_and_image(imgMoveY, input_info.dirO + cases[i].dir + "move/" + "y/" + cases[i].basename + ".raw");
		mhdRef.save_mhd_and_image(imgMoveZ, input_info.dirO + cases[i].dir + "move/" + "z/" + cases[i].basename + ".raw");

		//imgI�ɕ����Ǘ��ǂݍ����imgO��
		nari::vector<unsigned short> imgI(xeFl * yeFl * zeFl), imgO(xeRef * yeRef * zeRef);
		mhdFl.load_mhd_and_image(imgI, input_info.dirFl + cases[i].dir + cases[i].basename + ".mhd");
		//�����Ǘ����Ǘ�Ɉʒu���킹,�ό`��̌��摜��ۑ�
		deformation_using_movement(xeRef, yeRef, zeRef, xeFl, yeFl, zeFl, imgMoveX, imgMoveY, imgMoveZ, imgI, imgO);
		mhdRef.save_mhd_and_image(imgO, input_info.dirO + cases[i].dir + "org/" + cases[i].basename + ".raw");

		//���x����ό`���Ă���
		nari::vector<unsigned char> imgLabel;
		nari::mhd mhdLabel;
		std::ostringstream oss;
		oss << input_info.dirFl << cases[i].dir + cases[i].basename << "_label.mhd";

		std::cout <<oss.str() << std::endl;
		if (!nari::system::file_is_exist(oss.str())) continue;
		mhdLabel.load_mhd_and_image(imgLabel, oss.str());
		mhdLabel.size123(xeFl, yeFl, zeFl);

		// �����ϊ�
		nari::vector<float> imgLabelDist(xeFl * yeFl * zeFl), imgLabelDistO(xeRef * yeRef * zeRef);
		nari::distance::euclidean_signed_distance_transform(imgLabel.ptr(), imgLabelDist.ptr(), xeFl, yeFl, zeFl);
		mhdRef.save_mhd_and_image(imgLabelDist, input_info.dirO + cases[i].dir + "LS/" + cases[i].basename + ".raw");
		std::cout << "(�M�E�ցE�L)" << std::endl;
		//�ό`
		deformation_using_movement(xeRef, yeRef, zeRef, xeFl, yeFl, zeFl, imgMoveX, imgMoveY, imgMoveZ, imgLabelDist, imgLabelDistO);
		oss.str("");

		std::cout << "(�M�E�ցE�L)" << std::endl;

		//�����摜��臒l����
		nari::vector<unsigned char> imgLabelO(xeRef * yeRef * zeRef, 0);
		for (int s = 0; s < xeRef * yeRef * zeRef; s++)
		{
			if (imgLabelDistO[s] < 0) imgLabelO[s] = 1;
		}
		oss.str("");

		std::cout << "(�M�E�ցE�L)" << std::endl;

		//���x���̕ۑ�
		oss << input_info.dirO << cases[i].dir + "normalized_label/" + input_info.caseRef_name << "_nmlzd.raw";
		mhdRef.save_mhd_and_image(imgLabelO, oss.str());



	}
}