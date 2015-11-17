/*先頭で日本語を打ち込んでおけばソースツリーで表示したときに文字化けしないらしいので*/
#ifndef __TEMPLATE_MATCHONG__
#define __TEMPLATE_MATCHONG__

#include "narivectorpp.h"
#include "nmi_matching_It.h"
#include <omp.h>

template <class T , class M>
void template_mathcing ( const nari::vector<T> &imgRef, const nari::vector<M> &imgFl ,const nari::vector<nari::vector<int>> &DispRef,
	nari::vector<nari::vector<int>> &DispFl,int xeRef, int yeRef, int zeRef,int xeFl, int yeFl, int zeFl, 
	int tmp ,int rangex ,int rangey, int rangez)
{
	int disp_num = DispRef.size();
	
	int tmp_size = (tmp*2+1)*(tmp*2+1)*(tmp*2+1);

	//読み見込んだ計測点すべてを走査するループ
	#pragma omp parallel for schedule(dynamic) num_threads(2)
	for (int a = 0 ; a < disp_num ; a++ ){
		nari::vector<int> disp(3);
		nari::vector<unsigned short> tmp_Ref(tmp_size);
		int t = 0;
		int v = 0;
		//テンプレートマッチング前にあらかじめ位置合わせしたい側の画像のテンプレート作成
		for (int r = 0 ; r < 2*tmp+1 ; r++ ){
			for ( int q = 0 ; q < 2*tmp+1 ; q++ ){
				for (int p = 0 ; p < 2*tmp+1 ; p++ ){
					int x = DispRef[a][0] - tmp + p;
					int y = DispRef[a][1] - tmp + q;
					int z = DispRef[a][2] - tmp + r;
					//テンプレートが画像からはみ出た場合は折り返した画像を入れる
					if ( x > xeRef-1 ) x = 2*(xeRef-1) -x;
					if ( y > yeRef-1 ) y = 2*(yeRef-1) -y;
					if ( z > zeRef-1 ) z = 2*(zeRef-1) -z;
					int s = xeRef*yeRef*z + xeRef*y + x;
					//テンプレート内の画素を0～32に正規化
					//tmp_Ref[t] = imgRef[s]*32/65535;
					tmp_Ref[t] = imgRef[s];
					t++;

				}
			}
		}

		std::cout << "(^^)<テンプレート作った" <<std::endl;

		//評価値最大値と対応点の座標を格納する変数を定義

		//double nmi_max = 0;
		//相関係数を格納する変数を定義
		double cc = 0, cc_max = 0;
		int xs,ys,zs;

		//テンプレートマッチング開始
		for (int k = tmp; k < rangez * 2+1 -tmp ; k++){
			for (int j = tmp; j < rangey * 2+1 -tmp ; j++){
				for (int i = tmp; i < rangex * 2+1 -tmp ; i++){
					nari::vector<unsigned short> tmp_Fl(tmp_size);
					int u=0;
					//位置合わせされる側の画像テンプレート作成
					for (int r = 0 ; r < 2*tmp+1 ; r++ ){
						for (int q = 0 ; q < 2*tmp+1 ; q++ ){
							for (int p = 0 ; p < 2*tmp+1 ; p++ ){
								int x = DispRef[a][0] - rangex + i - tmp + p;
								int y = DispRef[a][1] - rangey + j - tmp + q;
								int z = DispRef[a][2] - rangez + k - tmp + r;
								//テンプレートが画像からはみ出た場合は折り返した画像を入れる
								if ( x > xeFl-1 ) x = 2*(xeFl-1) -x;
								if ( y > yeFl-1 ) y = 2*(yeFl-1) -y;
								if ( z > zeFl-1 ) z = 2*(zeFl-1) -z;
								int s = xeFl*yeFl*z + xeFl*y + x;
								//tmp_Fl[u] = imgFl[s]*32/65535;
								tmp_Fl[u] = imgFl[s];
								u++;

							}
						}
					}

					//平均をだす
					double meanref = 0.0, meanfl =0.0;
					for ( int c=0 ; c<tmp_size ; c++){
						meanref += tmp_Ref[c]/tmp_size;
						meanfl += tmp_Fl[c]/tmp_size;
					}

					double stdref = 0.0,stdfl = 0.0,cov = 0.0;
					for ( int c=0 ; c<tmp_size ; c++){
						stdref += (tmp_Ref[c]-meanref)*(tmp_Ref[c]-meanref);
						stdfl  += (tmp_Fl[c]-meanfl)*(tmp_Fl[c]-meanfl);
						cov    += (tmp_Ref[c]-meanref)*(tmp_Fl[c]-meanfl);
					}

					cc = abs(cov/(sqrt(stdref)*sqrt(stdfl)));

					if (cc > cc_max){
						cc_max = cc;
						xs = DispRef[a][0] - rangex + i;
						ys = DispRef[a][1] - rangey + j;
						zs = DispRef[a][2] - rangez + k;
					}

					////NMIを評価
					//double nmi = calc_NMI(tmp_Ref, tmp_Fl); //NMI（基準症例vs浮動症例）

					//完全に一致したときだけ
					/*int count = 0;
					for ( int c=0 ; c<u ; c++){
						int X = tmp_Ref[c]-tmp_Fl[c];
						if (X != 0)count++;
					}
					if (count == 0){
						xs = DispRef[a][0] - rangex + i;
						ys = DispRef[a][1] - rangey + j;
						zs = DispRef[a][2] - rangez + k;
					}*/

					
					//現在最大のnmiより大きければ対応点の座標を更新
					/*if (nmi > nmi_max){
						nmi_max = nmi;
						xs = DispRef[a][0] - rangex + i;
						ys = DispRef[a][1] - rangey + j;
						zs = DispRef[a][2] - rangez + k;
					}*/
					/*std::cout << "a=" << a <<std::endl;
					std::cout<<"i= "<<i<<" j= "<<j<<" k= "<<k<<std::endl; 
					system("pause");*/
				}

			}
			std::cout << "(´・ω・)<ここか？？" <<std::endl;
		}
		disp[0] = xs;
		disp[1] = ys;
		disp[2] = zs;
		DispFl.push_back(disp);

		std::cout << "(^^)<探索点みつけた" <<std::endl;
		std::cout << "a=" << a <<std::endl;
		std::cout << "x=" << DispRef[a][0] << " y=" << DispRef[a][1] <<" z=" << DispRef[a][2] <<std::endl;
		std::cout << "xs=" << xs << " ys=" << ys <<" zs=" << zs <<std::endl;
	}
}
#endif