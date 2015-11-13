#ifndef __TEMPLATE_MATCHONG__
#define __TEMPLATE_MATCHONG__

#include "narivectorpp.h"
#include "nmi_matching_It.h"


template <class T , class M>
void template_mathcing ( const nari::vector<T> &imgRef, const nari::vector<M> &imgFl ,const nari::vector<int> &x_disp ,
						const nari::vector<int> &y_disp ,const nari::vector<int> &z_disp ,nari::vector<int> &x_disp2 ,
						nari::vector<int> &y_disp2 ,nari::vector<int> &z_disp2 ,int xeRef, int yeRef, int zeRef,
						int xeFl, int yeFl, int zeFl, int tmp ,int rangex ,int rangey, int rangez)
{
	int x_num = x_disp.size();
	int y_num = y_disp.size();
	int z_num = z_disp.size();


	int tmp_size = (tmp*2+1)*(tmp*2+1)*(tmp*2+1);

	nari::vector<short> tmp_Ref(tmp_size);
	nari::vector<short> tmp_Fl(tmp_size);

	int a,b,c,p,q,r,s,x,y,z,t,u,v;


		//読み見込んだ計測点すべてを走査するループ
		for ( c = 0 ; c < z_num ; c++ ){
			for ( b = 0 ; b < y_num ; b++ ){
				for ( a = 0 ; a < x_num ; a++ ){
					t = 0;
					v = 0;
					//テンプレートマッチング前にあらかじめ位置合わせしたい側の画像のテンプレート作成
					for ( q = 0 ; q < 2*tmp+1 ; q++ ){
						for (  p = 0 ; p < 2*tmp+1 ; p++ ){
							for ( r = 0 ; r < 2*tmp+1 ; r++ ){
								x = x_disp[a] - tmp + r;
								y = y_disp[b] - tmp + q;
								z = z_disp[c] - tmp + p;
								//テンプレートが画像からはみ出た場合は折り返した画像を入れる
								if ( x > xeRef-1 ) x = 2*(xeRef-1) -x;
								if ( y > yeRef-1 ) y = 2*(yeRef-1) -y;
								if ( z > zeRef-1 ) z = 2*(zeRef-1) -z;
								s = xeRef*yeRef*z + xeRef*y + x;
								tmp_Ref[t] = imgRef[s];
								//テンプレート内の画素を0～32に正規化
								if( tmp_Ref[t] < 0) tmp_Ref[t] = 0;
								tmp_Ref[t] = tmp_Ref[t]*32 / 32767;
								t++;

							}
						}
					}

					//評価値最大値と対応点の座標を格納する変数を定義
					double nmi_max = 0;
					int xs,ys,zs;

					//テンプレートマッチング開始
					for (int k = tmp; k < rangez * 2 -tmp ; k++){
						for (int j = tmp; j < rangey * 2 -tmp ; j++){
							for (int i = tmp; i < rangex * 2 -tmp ; i++){
								u=0;
								//位置合わせされる側の画像テンプレート作成
								for ( r = 0 ; r < 2*tmp+1 ; r++ ){
									for ( q = 0 ; q < 2*tmp+1 ; q++ ){
										for ( p = 0 ; p < 2*tmp+1 ; p++ ){
											x = x_disp[a] - rangex + i - tmp + p;
											y = y_disp[b] - rangey + j - tmp + q;
											z = z_disp[c] - rangez + k - tmp + r;
											//テンプレートが画像からはみ出た場合は折り返した画像を入れる
											if ( x > xeFl-1 ) x = 2*(xeFl-1) -x;
											if ( y > yeFl-1 ) y = 2*(yeFl-1) -y;
											if ( z > zeFl-1 ) z = 2*(zeFl-1) -z;
											s = xeRef*yeRef*z + xeRef*y + x;
											tmp_Fl[u] = imgFl[s];
											//テンプレート内の画素を0～32に正規化
											if( tmp_Fl[u] < 0) tmp_Fl[u] = 0;
											tmp_Fl[u] = tmp_Fl[u]*32 / 32767;
											u++;

										}
									}
								}

								//NMIを評価
								double nmi = calc_NMI(tmp_Ref, tmp_Fl); //NMI（基準症例vs浮動症例）

								//現在最大のnmiより大きければ対応点の座標を更新
								if (nmi > nmi_max){
									nmi_max = nmi;
									xs = x_disp[a] - rangex + i;
									ys = y_disp[a] - rangey + j;
									zs = z_disp[a] - rangez + k;
								}
							}
						}
					}
					x_disp2[v] = xs;
					y_disp2[v] = ys;
					z_disp2[v] = zs;
					v++;
				}
			}
		}
}
#endif