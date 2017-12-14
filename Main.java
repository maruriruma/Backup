package b04;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Objects;

//基本的にコメントアウトで使いたいメカニズムを変更します。GUIは作ってません(><)
public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Main().start();
	}

	int conNum = 0;// 確定者の数計算
	final double siki = 0;
	double win;
	double envyin = 0; // 嫉妬人数(sum)
	double envyout = 0; // 嫉妬人数(sum)

	void start() {

		final int n = 512; // 学生数
		final int m = 64; // 学校数
		final int aR = 64; // 要素的下限
		final double a = 0.6; // alpha
		final double b = 0; // beta
		final int q = 40; // 学校個別の上限
		win=0;
		envyin = 0; // 嫉妬人数(sum)
		envyout = 0; // 嫉妬人数(sum)

		ArrayList<Student> Sts = new ArrayList<Student>(); // 学生の集合
		ArrayList<College> Cs = new ArrayList<College>(); // 学校の集合
		ArrayList<Region> Rs = new ArrayList<Region>(); // 地域の集合
		final int acq = 8; // ACDAの時のキャップ
		BigDecimal[] bdS = new BigDecimal[7];
		BigDecimal[] bdA = new BigDecimal[7];
		BigDecimal[] bdC = new BigDecimal[7];

		try {
			FileWriter fw = new FileWriter("C:/Users/admin/Desktop/Exceldata/test.csv", false);
			PrintWriter pw = new PrintWriter(new BufferedWriter(fw));

			for (int ct = 1; ct <= 7; ct++) {
				double envySum = 0; // 嫉妬人数(sum)
				double wasteSum = 0; // 空席要求人数(sum)
				win=0;



				final int cycle = 100; // シミュレーション回数
				double cSum1 = 0;
				double cSum2 = 0;
				int[] withinNumSum = new int[25]; // 累積分布関数に使用(sum)
				double[] averageRankSum = new double[m]; // 学校側からの学生のランクの平均(sum)
				double average = 0;
				ArrayList<Integer> qrs = new ArrayList<Integer>(); // 苦し紛れの策(
																	// ;∀;)
				ArrayList<Integer> prs = new ArrayList<Integer>(); // 苦し紛れの策(
																	// ;∀;)
				int[] complainS = new int[n];
				int maxC = 0;
				int minC = 600;
				int[] hozon = new int[cycle];

				for (int r = 0; r < cycle; r++) {
					// System.out.println((r + 1) + "回目");
					Arrays.fill(complainS, 0);
					// if (r == 0) {
					Sts = new ArrayList<Student>(); // 学生の集合
					Cs = new ArrayList<College>(); // 学校の集合
					Rs = new ArrayList<Region>(); // 地域の集合

					// setUp(n, m, a, b, Sts, Cs, q, aR, Rs);
					setUp(n, m, a, b, Sts, Cs, q, ct * aR, Rs);

//					ArrayList<Student> sA=new ArrayList<Student>();
//					int p=(int)Math.random()*512;
//					sA.add(sts.get(p));
//					for(Student s:Sts){
//						if(Math.random()<siki)
//							sA.add(s);
					//}
//					System.out.println("size: "+sA.size());
//					for (int i = 0; i < sA.size(); i++) {
//						Student s = sA.get(i);
/*						for (int j = 0; j < Cs.size(); j++) {
//							s.preS.get(j).pop += Cs.size() - j;
//						}
//					}
					ArrayList<College> clcs = new ArrayList<College>(Cs);
					PLmergeSort(clcs);
					int sum=0;
					for(Student s:Sts){
						for(College c:s.preS){
							sum+=Math.abs(clcs.indexOf(c)-s.preS.indexOf(c));
						}
					}
					System.out.println(sum);
					System.exit(0);

*/					Rs.get(0).qr = n;// Cの上限を学生数に
					// }
					qrs = new ArrayList<Integer>(); // 苦し紛れの策( ;∀;)
					prs = new ArrayList<Integer>(); // 苦し紛れの策( ;∀;)
					for (int i = 0; i < Rs.size(); i++) {
						qrs.add(Rs.get(i).qr);
						prs.add(Rs.get(i).pr);
					}
					//showStudentPre(Sts); //
					// 学生の好みが見たいときにコメントを外してください
					//showCollegePre(Cs); //
					// minS = "";
					for (int i = 0; i < Rs.size(); i++) {
						Rs.get(i).pr = prs.get(i);
						Rs.get(i).qr = qrs.get(i);
					}
					int[] withinNum = new int[25]; // 好みの何位以上に割り当てられたか(人数)
					int[] averageRank = new int[Cs.size()]; // 割り当てられた学生の合計優先順位

					// showStudentPre(Sts); //
					// 学生の好みが見たいときにコメントを外してください
					// showCollegePre(Cs); //
					// 学校の好みが見たいときにコメントを外してください
					// System.out.println("");

					long sta = System.currentTimeMillis();

					// 試したいメカニズムのコメントを外してください
					// DA(Sts, Cs, Rs, q);
					// ACDA(Sts, Cs, Rs, acq);
					// SDRQ(Sts, Cs, Rs, Sts);
					// MSDARQ(Sts, Cs, Rs,q);
					// ADA(Sts, Cs, Rs, q);
					// imADA(Sts, Cs, Rs,
					// q);//謝罪もの(MLの学生順を学校からの人気順にしたもの)
					// MACDA(Sts, Cs, Rs, q);
					//PLDARQ(Sts, Cs, Rs, Sts, Sts);
					// MSPLDARQ(Sts, Cs, Rs);
					multiPLDA(Sts, Cs, Rs);

					long end = System.currentTimeMillis();
					// System.out.println((end -
					// sta) + "ms");

					//showDetailMatch(Cs); //
					// どの学校にどの学生が割り当てられたか確認(詳細?)
					// showRoughMatch(Cs); //
					// どの学校にどの学生が割り当てられたか確認(おおざっぱ?)

					// 好みの何位以上に割り当てられたか計算
					for (int i = 0; i < withinNum.length; i++) {
						for (int j = 0; j < Sts.size(); j++) {
							Student s = Sts.get(j);
							if (s.preS.indexOf(s.c) <= i)
								withinNumSum[i]++;
						}
					}

					for (int i = 0; i < averageRank.length; i++) {
						if (Cs.get(i).hasS.size() == 0)
							continue;
						for (int j = 0; j < Cs.get(i).hasS.size(); j++) {
							if (Cs.get(i).hasS.size() == 0) {
								continue;
							}
							averageRank[i] += Cs.get(i).preC.indexOf(Cs.get(i).hasS.get(j)) + 1;
						}
						averageRankSum[i] += ((double) averageRank[i] / Cs.get(i).hasS.size());
					}

					envySum += envy(Sts, complainS); // 妥当な不満を持つ学生の人数を計算
					wasteSum += waste(qrs, prs, Sts, Rs, complainS); // 空きシートを要求する学生の人数を計算
					int counter = 0;
					for (int cs : complainS) {
						if (cs >= 1) {
							cSum1++;
							counter++;
						}
						if (cs == 2)
							cSum2++;
					}
					if (counter > maxC)
						maxC = counter;
					if (counter < minC)
						minC = counter;
					hozon[r] = counter;

				}
				double sigma = calcSigma(hozon, cSum1 / cycle);
				System.out.println("ok");
//				bdS[ct - 1] = new BigDecimal(envySum / cycle).setScale(3, BigDecimal.ROUND_HALF_UP);
//				bdA[ct - 1] = new BigDecimal(wasteSum / cycle).setScale(3, BigDecimal.ROUND_HALF_UP);
//				bdC[ct - 1] = new BigDecimal(cSum1 / cycle).setScale(3, BigDecimal.ROUND_HALF_UP);

				bdS[ct - 1] = new BigDecimal(envyin / cycle).setScale(3, BigDecimal.ROUND_HALF_UP);
				bdA[ct - 1] = new BigDecimal(envyout / cycle).setScale(3, BigDecimal.ROUND_HALF_UP);
				bdC[ct - 1] = new BigDecimal(win / cycle).setScale(3, BigDecimal.ROUND_HALF_UP);

				// BigDecimal bd=new BigDecimal(cSum1/cycle).setScale(3,
				// BigDecimal.ROUND_HALF_UP);
				// pw.print(bd.doubleValue());
				// pw.print(",");
				System.out.println("下限数:" + ct * aR);
				System.out.println("envy:" + envySum + "=" + (envySum / cycle));
				System.out.println("envyin:" + envyin + "=" + (envyin / cycle));
				System.out.println("envyout:" + envyout + "=" + (envyout / cycle));
				System.out.println("waste:" + wasteSum + "=" + (wasteSum / cycle));
				System.out.println("win:" + win + "=" + (win / cycle));
				System.out.println("conNum:" + conNum + "=" + ((double) conNum / cycle));
				System.out.println("重複考慮した場合の不満人数:" + cSum1 + "=" + (cSum1 / cycle));
				System.out.println("重複してる人数:" + cSum2 + "=" + (cSum2 / cycle));
				System.out.println("sigma:" + sigma);

				for (int i = 0; i < withinNumSum.length; i++) {
					System.out.println("within:" + (i + 1) + "=" + ((double) withinNumSum[i] / cycle));
				}

				for (int i = 0; i < averageRankSum.length; i++) {
					average += (averageRankSum[i] / cycle);
				}
				System.out.println("average = " + (average / m));
			}
			pw.print("envy" + ",");
			for (BigDecimal bd : bdS) {
				pw.print(bd.doubleValue() + ",");
			}
			pw.print("\nwaste" + ",");
			for (BigDecimal bd : bdA) {
				pw.print(bd.doubleValue() + ",");
			}
			pw.print("\nchou" + ",");
			for (BigDecimal bd : bdC) {
				pw.print(bd.doubleValue() + ",");
			}
			pw.close();
		} catch (IOException e) {
			// TODO 自動生成された catch ブロック
			e.printStackTrace();
		}
	}

	// セットアップ
	void setUp(int n, int m, double a, double b, ArrayList<Student> Sts, ArrayList<College> Cs, int q, int aR,
			ArrayList<Region> Rs) {
		Student.sum = 0;
		College.sum = 0;
		Student.m = m;
		College.n = n;
		int ReNum = m * 2 - 1; // 地域数
		Region.aR = aR; // 要素的下限制約数
		Region.aN = Region.aR / ReNum + 1; // 1地域辺りの要素的下限数
		Region.moreRe = Region.aR % ReNum;
		Region.lessRe = ReNum - Region.moreRe;

		for (int i = 0; i < n; i++)
			Sts.add(new Student(Math.random())); // 学生生成
		for (int i = 0; i < m; i++)
			Cs.add(new College(Math.random())); // 学校生成
		for (int i = 0; i < n; i++)
			Sts.get(i).setUp(a, Cs); // 学生好み生成
		for (int i = 0; i < m; i++)
			Cs.get(i).setUp(b, Sts); // 学校好み生成

		ArrayList<Region> stack = new ArrayList<Region>();
		stack.add(new Region(null, Cs, q, 0));

		while (!stack.isEmpty()) { // 階層的地域作成
			Region node = stack.get(0);
			stack.remove(0);
			Rs.add(node);
			int cycle = 0;
			for (Region childNode : node.childrenN) {
				stack.add(cycle, childNode);
				cycle++;
			}
		}
		//showRegionNum(Rs); // 地域の詳細(コメント外してください)
	}

	// DA(返り値をマッチングとかにした方が良いのかなー)
	void DA(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs, int aq) {

		ArrayList<Integer> fixNum = new ArrayList<Integer>(); // MSDAQRで使う(固定した前までのマッチングを保存する)

		for (int i = 0; i < cs.size(); i++) { // 各学校の保存生徒数を格納
			fixNum.add(cs.get(i).hasS.size());
		}

		ArrayList<College> Xc = new ArrayList<College>();
		ArrayList<Integer> Xnum = new ArrayList<Integer>();

		while (!(isAllStuMatch(sts))) {// DA
			for (int i = 0; i < cs.size(); i++) {// 固定されている学生以外を削除
				int max = cs.get(i).hasS.size();
				for (int j = 0; j < max - fixNum.get(i); j++) {
					cs.get(i).hasS.remove(cs.get(i).hasS.size() - 1);
				}
			}

			// 木のアップデートで使う
			Xc = new ArrayList<College>();
			Xnum = new ArrayList<Integer>();

			for (int i = 0; i < sts.size(); i++) {// 学生が学校に応募
				Student s = sts.get(i);
				while (true) {
					if (s.preS.get(s.maxPre).forbidden)
						s.maxPre++;
					else
						break;
				}
				s.preS.get(s.maxPre).hasS.add(s);
				s.c = s.preS.get(s.maxPre);
				if (!Xc.contains(s.c)) {
					Xc.add(s.c);
					Xnum.add(1);
				} else {
					Xnum.set(Xc.indexOf(s.c), Xnum.get(Xc.indexOf(s.c)) + 1);
				}
			}

			for (int i = 0; i < cs.size(); i++) {// 好み順に並び替えて上限以上の学生を拒否
				College c = cs.get(i);
				int fix = fixNum.get(i);
				mergeSort(fix, c.hasS, c);
				if (c.hasS.size() > aq) {
					int rem = c.hasS.size() - aq;
					for (int j = 0; j < rem; j++) {
						c.hasS.get(aq).c = null;
						c.hasS.get(aq).maxPre++;
						c.hasS.remove(aq);
					}
				}
			}
		}
		update(Xc, Xnum, rs);

	}

	void ACDA(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs, int acq) {
		DA(sts, cs, rs, acq);
	}

	void SDRQ(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs, ArrayList<Student> ML) {

		for (int k = 0; k < sts.size(); k++) {
			int pC = 0;
			Student s = ML.get(k);
			boolean already = false;

			ArrayList<College> Xc = new ArrayList<College>();
			ArrayList<Integer> Xnum = new ArrayList<Integer>();

			College c = new College(0);

			for (int i = 0; i < rs.size(); i++) {// 根ノードの下限計算
				pC += rs.get(i).ar;
			}

			if (pC < sts.size() - k) {// 根ノードの下限が残り人数より多い
				for (int i = 0; i < cs.size(); i++) {
					c = s.preS.get(s.maxPre);

					for (int j = 0; j < rs.size() && !already; j++) {
						if (rs.get(j).hasC.size() > 1)
							continue;
						if (rs.get(j).hasC.get(0).equals(c)) {
							if (rs.get(j).qr > 0) {
								s.preS.get(s.maxPre).hasS.add(s);
								s.c = s.preS.get(s.maxPre);
								already = true;
							}
						}
					}
					if (already)
						break;
					s.maxPre++;
				}
			} else if (pC == sts.size() - k) {// 根ノードの下限が残り人数と一緒
				for (int i = 0; i < cs.size(); i++) {
					c = s.preS.get(s.maxPre);

					for (int j = 0; j < rs.size() && !already; j++) {
						if (rs.get(j).hasC.size() > 1)
							continue;
						if (rs.get(j).hasC.get(0).equals(c)) {
							if (rs.get(j).qr > 0 && Ear(c, rs)) {
								s.preS.get(s.maxPre).hasS.add(s);
								s.c = s.preS.get(s.maxPre);
								already = true;
							}
						}
					}
					if (already)
						break;
					s.maxPre++;
				}
			}
			Xc.add(c);
			Xnum.add(1);
			update(Xc, Xnum, rs);
		}
	}

	void MSDARQ(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs, int q) {

		ArrayList<Student> ML = (ArrayList<Student>) sts.clone();
		ArrayList<Student> E = (ArrayList<Student>) ML.clone();

		int count = 0;
		for (int k = 0; k < sts.size() && !(isAllStuMatch(sts)); k++) {
			int e = rs.get(0).pr;
			count++;

			ArrayList<Student> Ek = new ArrayList<Student>();
			for (int i = 0; i < e; i++) {
				Ek.add(E.get(E.size() - e + i));
			}
			ArrayList<Student> ml = new ArrayList();

			E.removeAll(Ek);

			if (E.size() != 0) {
				DA(E, cs, rs, q);
			} else {
				SDRQ(Ek, cs, rs, ml);
			}
			E = Ek;
		}
	}

	void ADA(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs, int q) {

		ArrayList<Student> L = new ArrayList<Student>(sts);
		ArrayList<Integer> prs = new ArrayList<Integer>();// 下限保存
		ArrayList<Integer> qrs = new ArrayList<Integer>();// 上限保存

		for (int i = 0; i < rs.size(); i++) {
			prs.add(rs.get(i).pr);
			qrs.add(rs.get(i).qr);
		}
		while (L.size() > 0 && !isAllStuMatch(sts)) {// 残りの学生がいる、かつまだ全員が割り当てられてない場合続く
			int t = 0;
			ArrayList<Student> DAStu = new ArrayList<Student>();

			ArrayList<Integer> fixNum = new ArrayList<Integer>(); // MSDAQRで使う(固定した前までのマッチングを保存する)

			for (int i = 0; i < cs.size(); i++) { // 各学校の保存生徒数を格納
				fixNum.add(cs.get(i).hasS.size());
			}

			while (t < L.size()) {
				ArrayList<Region> subRe = new ArrayList<Region>();

				for (int i = 0; i < cs.size(); i++) {
					int max = cs.get(i).hasS.size();
					for (int j = 0; j < max - fixNum.get(i); j++) {
						cs.get(i).hasS.remove(cs.get(i).hasS.size() - 1);
					}
				}

				DAStu = new ArrayList<Student>();

				ArrayList<Region> stack = new ArrayList<Region>();
				stack.add(new Region(null, cs, q, 0));
				while (!stack.isEmpty()) { // 階層的地域作成
					Region node = stack.get(0);
					stack.remove(0);
					subRe.add(node);
					int cycle = 0;
					for (Region childNode : node.childrenN) {
						stack.add(cycle, childNode);
						cycle++;
					}
				}
				for (int i = 0; i < subRe.size(); i++) {
					subRe.get(i).pr = rs.get(i).pr;
					subRe.get(i).qr = rs.get(i).qr;
					subRe.get(i).ar = rs.get(i).ar;
				}

				for (int i = 0; i < t + 1; i++) {// Lのt番目の学生までDAする学生のリストに格納
					DAStu.add(L.get(i));
				}
				DA(DAStu, cs, subRe, q);

				if (subRe.get(0).pr < L.size() - (t + 1)) {// 制約に違反しないなら(階層的地域なので下限のみ参照)
					t++;
					continue;
				} else if (subRe.get(0).pr == L.size() - (t + 1)) {// 制約に違反するなら
					break;
				} else {
					System.exit(1);
					return;
				}
			}
			for (int i = 0; i < cs.size(); i++) {
				int max = cs.get(i).hasS.size();
				for (int j = 0; j < max - fixNum.get(i); j++) {
					cs.get(i).hasS.remove(cs.get(i).hasS.size() - 1);
				}
			}
			for (int i = 0; i < DAStu.size(); i++) {
				DAStu.get(i).c = null;
			}
			DA(DAStu, cs, rs, q);
			for (int i = 0; i < cs.size(); i++) {
				if (cs.get(i).forbidden)
					continue;
				else if (!Ear(cs.get(i), rs)) // その学校にいれると制約に違反してしまう
					cs.get(i).forbidden = true;
			}
			L.removeAll(DAStu);
		}
	}

	// 謝罪ものなので、このメソッドは気にしないでください
	void imADA(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs, int q) {
		ArrayList<Student> L = new ArrayList<Student>(sts);// 初期状態
		for (int i = 0; i < cs.size(); i++) {
			for (int j = 0; j < sts.size(); j++) {
				// cs.get(i).preC.get(j).pop += j;// 人気の平均
				if (cs.get(i).preC.get(j).pop == 0 || cs.get(i).preC.get(j).pop > j + 1) {
					cs.get(i).preC.get(j).pop = j + 1;// min人気
				}
			}
		}
		popmergeSort(L);// popular sort
		ArrayList<Integer> prs = new ArrayList<Integer>();
		ArrayList<Integer> qrs = new ArrayList<Integer>();
		for (int i = 0; i < rs.size(); i++) {
			prs.add(rs.get(i).pr);
			qrs.add(rs.get(i).qr);
		}
		while (L.size() > 0 && !isAllStuMatch(sts)) {
			int t = 0;
			ArrayList<Student> DAStu = new ArrayList<Student>();

			ArrayList<Integer> fixNum = new ArrayList<Integer>(); // MSDAQRで使う(固定した前までのマッチングを保存する)
			for (int i = 0; i < cs.size(); i++) { // 各学校の保存生徒数を格納
				fixNum.add(cs.get(i).hasS.size());
			}
			ArrayList<College> forCol = new ArrayList<College>();// forbbiden
																	// Colege

			while (t < L.size()) {
				for (int i = 0; i < cs.size(); i++) {
					int max = cs.get(i).hasS.size();
					for (int j = 0; j < max - fixNum.get(i); j++) {
						cs.get(i).hasS.remove(cs.get(i).hasS.size() - 1);
					}
				}

				DAStu = new ArrayList<Student>();
				for (int i = 0; i < t + 1; i++) {
					DAStu.add(L.get(i));
				}
				imDA(DAStu, cs, rs, q);
				int count = 0;
				for (int i = 0; i < cs.size(); i++) {
					College c = cs.get(i);
					if (c.forbidden)
						continue;
					if (c.hasS.size() >= q)
						continue;
					c.hasS.add(new Student());
					minimumCount(rs.get(0));
					for (int j = 0; j < rs.size(); j++) {
						if (rs.get(j).er > rs.get(j).qr) {
							count++;
							forCol.add(c);
							break;
						}
					}
					c.hasS.remove(c.hasS.size() - 1);
				}

				if (count == 0) {
					t++;
					continue;
				} else if (count > 0) {
					break;
				} else {
					System.exit(1);
					return;
				}
			}
			for (int i = 0; i < cs.size(); i++) {
				int max = cs.get(i).hasS.size();
				for (int j = 0; j < max - fixNum.get(i); j++) {
					cs.get(i).hasS.remove(cs.get(i).hasS.size() - 1);
				}
			}
			for (int i = 0; i < DAStu.size(); i++) {
				DAStu.get(i).c = null;
			}
			imDA(DAStu, cs, rs, q);
			for (int i = 0; i < forCol.size(); i++) {
				forCol.get(i).forbidden = true;
			}
			L.removeAll(DAStu);
		}
	}

	// 戦略的に操作可能だったもの、すごい時間かかる
	void MACDA(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs, int q) {

		final ArrayList<College> MLC = new ArrayList<College>(cs);
		sts.get(0).mergeSort(MLC);

		for (int i = 0; i < sts.size(); i++) {
			System.out.println("cycle" + i);
			ArrayList<Student> L = new ArrayList<Student>(sts);
			ArrayList<Integer> qrs = new ArrayList<Integer>();// 苦し紛れの策( ;∀;)
			ArrayList<Integer> prs = new ArrayList<Integer>();// 苦し紛れの策( ;∀;)
			ArrayList<Integer> ars = new ArrayList<Integer>();// 苦し紛れの策( ;∀;)

			for (int j = 0; j < rs.size(); j++) {
				qrs.add(rs.get(j).qr);
				prs.add(rs.get(j).pr);
				ars.add(rs.get(j).ar);
			}

			Student s = L.remove(i);
			L.add(s);

			ADA(L, cs, rs, q);

			for (int j = 0; j < cs.size(); j++) {
				int min = 0;
				if (!Objects.equals(cs.get(j).hasS, null))
					min = cs.get(j).hasS.size();
				if (cs.get(j).andNum > min)
					cs.get(j).andNum = min;
			}

			for (int j = 0; j < sts.size(); j++) {
				sts.get(j).c = null;
				sts.get(j).maxPre = 0;
			}

			for (int j = 0; j < cs.size(); j++) {
				cs.get(j).hasS = new ArrayList<Student>();
				cs.get(j).forbidden = false;
			}

			for (int j = 0; j < rs.size(); j++) {
				rs.get(j).pr = prs.get(j);
				rs.get(j).qr = qrs.get(j);
				rs.get(j).ar = ars.get(j);
			}
		}
		int sum = 0;
		int rem = 0;

		for (int i = 0; i < cs.size(); i++)
			sum += cs.get(i).andNum;
		if (sum < sts.size())
			rem = sts.size() - sum;
		System.out.println("rem:" + rem + ", sum:" + sum);
		int count = 0;
		while (rem > 0) {
			if (cs.get(count).andNum >= q) {
				count++;
				continue;
			}
			cs.get(count).andNum++;
			rem--;
		}
		sum = 0;
		for (int i = 0; i < cs.size(); i++)
			sum += cs.get(i).andNum;
		if (sum < sts.size())
			rem = sts.size() - sum;
		System.out.println("rem:" + rem + ", sum:" + sum);
		for (int i = 0; i < cs.size(); i++)
			System.out.println(i + "=" + cs.get(i).andNum);
		int[] andNum = new int[cs.size()];
		for (int i = 0; i < cs.size(); i++) {
			andNum[i] = cs.get(i).andNum;
		}
		DA(sts, cs, rs, andNum);
	}

	// MACDA用のDAだった(少し引数が異なるだけだと思います)
	void DA(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs, int[] andNum) {

		ArrayList<Integer> fixNum = new ArrayList<Integer>(); // MSDAQRで使う(固定した前までのマッチングを保存する)
		for (int i = 0; i < cs.size(); i++) { // 各学校の保存生徒数を格納
			fixNum.add(cs.get(i).hasS.size());
		}
		ArrayList<College> Xc = new ArrayList<College>();
		ArrayList<Integer> Xnum = new ArrayList<Integer>();
		while (!(isAllStuMatch(sts))) {// DA
			for (int i = 0; i < cs.size(); i++) {
				int max = cs.get(i).hasS.size();
				for (int j = 0; j < max - fixNum.get(i); j++) {
					cs.get(i).hasS.remove(cs.get(i).hasS.size() - 1);
				}
			}
			Xc = new ArrayList<College>();
			Xnum = new ArrayList<Integer>();

			for (int i = 0; i < sts.size(); i++) {
				Student s = sts.get(i);
				while (true) {
					if (s.preS.get(s.maxPre).forbidden)
						s.maxPre++;
					else
						break;
				}
				s.preS.get(s.maxPre).hasS.add(s);
				s.c = s.preS.get(s.maxPre);
				if (!Xc.contains(s.c)) {
					Xc.add(s.c);
					Xnum.add(1);
				} else {
					Xnum.set(Xc.indexOf(s.c), Xnum.get(Xc.indexOf(s.c)) + 1);
				}
			}
			for (int i = 0; i < cs.size(); i++) {
				College c = cs.get(i);
				int fix = fixNum.get(i);
				mergeSort(fix, c.hasS, c);
				if (c.hasS.size() > andNum[i]) {
					int rem = c.hasS.size() - andNum[i];
					for (int j = 0; j < rem; j++) {
						c.hasS.get(andNum[i]).c = null;
						c.hasS.get(andNum[i]).maxPre++;
						c.hasS.remove(andNum[i]);
					}
				}
			}
		}
		update(Xc, Xnum, rs);
	}

	// PLの順序を変えるときはメソッド内の最初の所でコメントでいじってください
	void PLDARQ(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs, ArrayList<Student> sA,
			ArrayList<Student> sB) {
		ArrayList<Contract> PL = new ArrayList<Contract>();

		// この下から使いたいPL順序を1つだけコメントを外してください。

		//tiebreakPLSet(PL, sA, cs); // タイブレーク順序
		// verticalTiebreakPLSet(PL, sts, cs); // 縦のタイブレークセット
		// preferFixPLSet(PL, sts, cs, rs); // 確定者第一志望のみ
		// preferAllPLSet(PL, sts, cs, rs); // 確定者全選好反映
		// preferViPLSet(PL, sts, cs, rs); // 確定者全選好反映(確定者が悪人)
		allStuPLSet(PL, sB, cs, rs, sA); // 全生徒選好反映
		//preStuPLSet(PL, sB, cs, rs, sA);

		// for(Contract pl:PL)
		// System.out.println(pl.getStudent().number+","+pl.getCollege().number);

		// ここまで

		minimumCount(rs.get(0));// 期待最小カウント

		ArrayList<Contract> X = new ArrayList<Contract>();// 契約の集合

		ArrayList<Student> rejectS = new ArrayList<Student>(sts);// 拒否された学生の集合

		while (!(isAllStuMatch(sts))) {
			int sum = 0;
			for (int i = 0; i < cs.size(); i++) {
				sum += cs.get(i).hasS.size();
			}

			// 初期化
			for (int i = 0; i < sts.size(); i++) {
				sts.get(i).c = null;
			}
			for (int i = 0; i < cs.size(); i++) {
				cs.get(i).hasS = new ArrayList<Student>();
			}

			for (int i = 0; i < rejectS.size(); i++) {
				Student s = rejectS.get(i);
				Contract x = new Contract(s, s.preS.get(s.maxPre));
				boolean fin = true;
				while (fin) {
					for (int j = 0; j < PL.size(); j++) {
						if (x.contractContain(PL.get(j))) {
							x.setPriority(PL.get(j).getPriority());
							fin = false;
							break;
						}
					}
				}
				X.add(x);
			}
			mergeSort(X);

			rejectS = new ArrayList<Student>();
			ArrayList<Contract> rejX = new ArrayList<Contract>();

			for (int i = 0; i < X.size(); i++) {
				Contract x = X.get(i);
				Student s = x.getStudent();
				College c = x.getCollege();
				c.hasS.add(s);
				minimumCount(rs.get(0));
				boolean ban = false;
				for (int j = 0; j < rs.size(); j++) {
					Region r = rs.get(j);
					if (r.er > r.qr) {// 学生sを学校cに追加したときに違反が発生
						c.hasS.remove(s);
						s.maxPre++;
						rejectS.add(s);
						rejX.add(x);
						ban = true;
						break;
					}
				}
				if (ban)
					continue;
				s.c = c;
			}
			X.removeAll(rejX);
		}
	}

	// これも謝罪ものだったので気にしないでください
	void MSPLDARQ(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs) {
		for (int i = 0; i < cs.size(); i++) {
			for (int j = 0; j < sts.size(); j++) {
				if (cs.get(i).preC.get(j).pop == 0 || cs.get(i).preC.get(j).pop > j + 1) {
					cs.get(i).preC.get(j).pop = j + 1;// min人気
				}
			}
		}
		ArrayList<Student> clst = new ArrayList<Student>(sts);
		popmergeSort(clst);
		minimumCount(rs.get(0));
		ArrayList<Student> greatS = new ArrayList<Student>();
		for (int i = 0; i < sts.size() - rs.get(0).pr; i++) {
			greatS.add(clst.get(0));
			clst.remove(0);
		}
		ArrayList<Contract> PL = new ArrayList<Contract>();
		tiebreakPLSet(PL, sts, cs);
		ArrayList<Contract> X = new ArrayList<Contract>();// 契約の集合
		ArrayList<Student> rejectS = new ArrayList<Student>(greatS);
		while (!(isAllStuMatch(greatS))) {
			int sum = 0;
			for (int i = 0; i < greatS.size(); i++) {
				greatS.get(i).c = null;
			}
			for (int i = 0; i < cs.size(); i++) {
				cs.get(i).hasS = new ArrayList<Student>();
			}
			for (int i = 0; i < rejectS.size(); i++) {
				Student s = rejectS.get(i);
				Contract x = new Contract(s, s.preS.get(s.maxPre));
				boolean fin = true;
				while (fin) {
					for (int j = 0; j < PL.size(); j++) {
						if (x.contractContain(PL.get(j))) {
							x.setPriority(PL.get(j).getPriority());
							fin = false;
							break;
						}
					}
				}
				X.add(x);
			}
			mergeSort(X);
			rejectS = new ArrayList<Student>();
			ArrayList<Contract> rejX = new ArrayList<Contract>();
			for (int i = 0; i < X.size(); i++) {
				Contract x = X.get(i);
				Student s = x.getStudent();
				College c = x.getCollege();
				c.hasS.add(s);
				minimumCount(rs.get(0));
				boolean ban = false;
				for (int j = 0; j < rs.size(); j++) {
					Region r = rs.get(j);
					if (r.er > r.qr) {
						c.hasS.remove(s);
						s.maxPre++;
						rejectS.add(s);
						rejX.add(x);
						ban = true;
						break;
					}
				}
				if (ban)
					continue;
				s.c = c;
			}
			X.removeAll(rejX);
		}
		ArrayList<College> clcs = new ArrayList<College>(cs);
		for (int i = 0; i < clcs.size(); i++) {
			clcs.get(i).pop = clcs.get(i).hasS.size();
		}
		PLmergeSort(clcs);
		minimumCount(rs.get(0));
		PL = new ArrayList<Contract>();
		greatSPLSet(PL, sts, clcs);
		X = new ArrayList<Contract>();// 契約の集合
		rejectS = new ArrayList<Student>(clst);
		while (!(isAllStuMatch(clst))) {
			int sum = 0;
			for (int i = 0; i < clst.size(); i++) {
				clst.get(i).c = null;
			}
			for (int i = 0; i < cs.size(); i++) {
				cs.get(i).hasS = new ArrayList<Student>();
			}
			for (int i = 0; i < rejectS.size(); i++) {
				Student s = rejectS.get(i);
				Contract x = new Contract(s, s.preS.get(s.maxPre));
				boolean fin = true;
				while (fin) {
					for (int j = 0; j < PL.size(); j++) {
						if (x.contractContain(PL.get(j))) {
							x.setPriority(PL.get(j).getPriority());
							fin = false;
							break;
						}
					}
				}
				X.add(x);
			}
			mergeSort(X);
			rejectS = new ArrayList<Student>();
			for (int i = 0; i < greatS.size(); i++) {
				Student s = greatS.get(i);
				s.c.hasS.add(s);
			}
			ArrayList<Contract> rejX = new ArrayList<Contract>();
			for (int i = 0; i < X.size(); i++) {
				Contract x = X.get(i);
				Student s = x.getStudent();
				College c = x.getCollege();
				c.hasS.add(s);
				minimumCount(rs.get(0));
				boolean ban = false;
				for (int j = 0; j < rs.size(); j++) {
					Region r = rs.get(j);
					if (r.er > r.qr) {
						c.hasS.remove(s);
						s.maxPre++;
						rejectS.add(s);
						rejX.add(x);
						ban = true;
						break;
					}
				}
				if (ban)
					continue;
				s.c = c;
			}
			X.removeAll(rejX);
		}
	}

	void multiPLDA(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs) {
		ArrayList<Student> sA = new ArrayList<Student>();
		ArrayList<Student> sB = new ArrayList<Student>();
		int b = sts.size() / 2;
		int a = sts.size() - b;

		int[] complainS = new int[sts.size()];

		for (int i = 0; i < sts.size(); i++) {
			if (a == 0)
				sB.add(sts.get(i));
			else if (b == 0)
				sA.add(sts.get(i));
			else {
				double r = Math.random();
				if (r < 0.5) {
					sA.add(sts.get(i));
					sts.get(i).group=0;
					a--;
				} else {
					sB.add(sts.get(i));
					sts.get(i).group=1;
					b--;
				}
			}
		}

/*		System.out.println("sA");
		for (Student s : sA) {
			System.out.println(s.number);
		}
		System.out.println("sB");
		for (Student s : sB) {
			System.out.println(s.number);
		}
*/
		ArrayList<ArrayList<Student>> ss = new ArrayList<ArrayList<Student>>();
		ArrayList<Integer> qrs = new ArrayList<Integer>();
		ArrayList<Integer> prs = new ArrayList<Integer>();
		ArrayList<Integer> qrA = new ArrayList<Integer>();
		ArrayList<Integer> prA = new ArrayList<Integer>();
		ArrayList<Integer> qrB = new ArrayList<Integer>();
		ArrayList<Integer> prB = new ArrayList<Integer>();
		for (int i = 0; i < rs.size(); i++) {
			qrs.add(rs.get(i).qr);
			prs.add(rs.get(i).pr);
			qrA.add(rs.get(i).qr / 2);
			prA.add(rs.get(i).pr / 2);
			rs.get(i).qr = qrA.get(i);
			rs.get(i).pr = prA.get(i);
		}
		int[] pS = new int[10];
		Arrays.fill(pS, 1000);
		for (Region r : rs) {
			if (pS[r.nest] == 1000)
				pS[r.nest] = r.pr;
			else
				pS[r.nest] += r.pr;
		}
//		for (int p : pS)
	//		System.out.println(p);

//		 for(Region r:rs){ System.out.println("rs"); System.out.println(r.qr);
//		 System.out.println(r.pr); }
		 PLDARQ(sA, cs, rs, sts, sB);
//		showDetailMatch(cs);
//		System.out.println("view");
		win+=waste2(qrA, prA, sA, rs);
		envyin+=envy2(sA);
		int count = 0;
		for (College c : cs) {
			ArrayList<Student> ns = new ArrayList<Student>();
			for (Student s : c.hasS) {
				ns.add(s);
			}
			ss.add(ns);
			count++;
			c.hasS.clear();
		}
		for (int i = 0; i < rs.size(); i++) {
			rs.get(i).qr = qrs.get(i) - qrA.get(i);
			rs.get(i).pr = prs.get(i) - prA.get(i);
			qrB.add(qrs.get(i) - qrA.get(i));
			prB.add(prs.get(i) - prA.get(i));
		}
		Arrays.fill(pS, 1000);
		for (Region r : rs) {
			if (pS[r.nest] == 1000)
				pS[r.nest] = r.pr;
			else
				pS[r.nest] += r.pr;
		}
//		for (int p : pS)
	//		System.out.println(p);

//		for(Region r:rs){ System.out.println("rs"); System.out.println(r.qr);
	//	System.out.println(r.pr); }
		PLDARQ(sB, cs, rs, sts, sA);
//		int wa=waste2(qrB, prB, sts, rs);
		win+=waste2(qrB, prB, sB, rs);
		envyin+=envy2(sB);
//		showDetailMatch(cs);
		for (int i = 0; i < cs.size(); i++) {
			for (Student s : ss.get(i)) {
				cs.get(i).hasS.add(s);
			}
		}
		for (int i = 0; i < rs.size(); i++) {
			rs.get(i).qr = qrs.get(i);
			rs.get(i).pr = prs.get(i);
		}

	}

	// タイブレーク順
	void tiebreakPLSet(ArrayList<Contract> PL, ArrayList<Student> sts, ArrayList<College> cs) {
		for (int i = 0; i < sts.size(); i++) {
			for (int j = 0; j < cs.size(); j++) {
				College c = cs.get(j);
				PL.add(new Contract(c.preC.get(i), c));
			}
		}
	}

	// 学校が全部の契約を格納してから次の学校に
	void verticalTiebreakPLSet(ArrayList<Contract> PL, ArrayList<Student> sts, ArrayList<College> cs) {
		for (int i = 0; i < cs.size(); i++) {
			College c = cs.get(i);
			for (int j = 0; j < sts.size(); j++) {
				PL.add(new Contract(c.preC.get(j), c));
			}
		}
	}

	// 確定者の第一志望のみ反映
	void preferFixPLSet(ArrayList<Contract> PL, ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs) {
		for (int i = 0; i < cs.size(); i++) {
			College c = cs.get(i);
			Region r = null;
			for (int j = 0; j < rs.size(); j++) {
				if (rs.get(j).hasC.size() != 1)
					continue;
				if (rs.get(j).hasC.get(0).equals(c)) {
					r = rs.get(j);
					break;
				}
			}
			for (int j = 0; j < r.pr; j++) {
				if (c.preC.get(j).preS.get(0).equals(c)) {
					c.pop++;
					conNum++;
				}
			}
		}
		System.out.println("確定者：" + conNum);
		ArrayList<College> clcs = new ArrayList<College>(cs);
		PLmergeSort(clcs);
		for (int i = 0; i < clcs.size(); i++) {
			College c = clcs.get(i);
			if (c.pop > 0) {
				for (int j = 0; j < sts.size(); j++) {
					PL.add(new Contract(c.preC.get(j), c));
				}
			} else
				break;
		}
		for (int i = 0; i < sts.size(); i++) {
			for (int j = 0; j < clcs.size(); j++) {
				College c = clcs.get(j);
				if (c.pop > 0) {
					continue;
				} else
					PL.add(new Contract(c.preC.get(i), c));
			}
		}
	}

	// 確定者の全選好反映
	void preferAllPLSet(ArrayList<Contract> PL, ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs) {
		ArrayList<Student> conS = new ArrayList<Student>();
		for (int i = 0; i < cs.size(); i++) {
			College c = cs.get(i);
			Region r = null;
			for (int j = 0; j < rs.size(); j++) {
				if (rs.get(j).hasC.size() != 1)
					continue;
				if (rs.get(j).hasC.get(0).equals(c)) {
					r = rs.get(j);
					break;
				}
			}
			for (int j = 0; j < r.pr; j++) {
				if (c.preC.get(j).preS.get(0).equals(c)) {
					conNum++;
					conS.add(c.preC.get(j));
				}
			}
		}
		System.out.println("確定者：" + conNum);
		if (conS.size() == 0) {
			tiebreakPLSet(PL, sts, cs);
			return;
		}
		System.out.println("確定者は0じゃない！");
		for (int i = 0; i < conS.size(); i++) {
			Student s = conS.get(i);
			for (int j = 0; j < cs.size(); j++) {
				s.preS.get(j).pop += cs.size() - j;
			}
		}

		ArrayList<College> clcs = new ArrayList<College>(cs);
		PLmergeSort(clcs);
		for (int i = 0; i < clcs.size(); i++) {
			College c = clcs.get(i);
			for (int j = 0; j < sts.size(); j++) {
				PL.add(new Contract(c.preC.get(j), c));
			}
		}
	}

	// 確定者の全選好反映(しかし確定者の選好を社会の学校の人気順の逆にする)
	void preferViPLSet(ArrayList<Contract> PL, ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs) {
		ArrayList<Student> conS = new ArrayList<Student>();
		for (int i = 0; i < cs.size(); i++) {
			College c = cs.get(i);
			Region r = null;
			for (int j = 0; j < rs.size(); j++) {
				if (rs.get(j).hasC.size() != 1)
					continue;
				if (rs.get(j).hasC.get(0).equals(c)) {
					r = rs.get(j);
					break;
				}
			}
			for (int j = 0; j < r.pr; j++) {
				if (c.preC.get(j).preS.get(0).equals(c)) {
					conNum++;
					conS.add(c.preC.get(j));
				}
			}
		}
		System.out.println("確定者：" + conNum);
		if (conS.size() == 0) {
			tiebreakPLSet(PL, sts, cs);
			return;
		}
		System.out.println("確定者は0じゃない！");
		for (int i = 0; i < sts.size(); i++) {
			Student s = sts.get(i);
			for (int j = 0; j < cs.size(); j++) {
				s.preS.get(j).pop += j;
			}
		}
		ArrayList<College> clcs = new ArrayList<College>(cs);
		PLmergeSort(clcs);
		for (int i = 0; i < conS.size(); i++) {
			ArrayList<College> copy = new ArrayList<College>(clcs);
			copy.remove(conS.get(i).preS.get(0));
			conS.get(i).preS.removeAll(copy);
			conS.get(i).preS.addAll(copy);
		}
		for (int i = 0; i < cs.size(); i++) {
			cs.get(i).pop = 0;
		}
		for (int i = 0; i < conS.size(); i++) {
			Student s = conS.get(i);
			for (int j = 0; j < cs.size(); j++) {
				s.preS.get(j).pop += cs.size() - j;
			}
		}
		PLmergeSort(clcs);
		for (int i = 0; i < clcs.size(); i++) {
			College c = clcs.get(i);
			for (int j = 0; j < sts.size(); j++) {
				PL.add(new Contract(c.preC.get(j), c));
			}
		}
	}

	// 全学生の全選好反映、もちろん戦略的操作可能
	void allStuPLSet(ArrayList<Contract> PL, ArrayList<Student> sA, ArrayList<College> cs, ArrayList<Region> rs,
			ArrayList<Student> sts) {
		for (int i = 0; i < sA.size(); i++) {
			Student s = sA.get(i);
			for (int j = 0; j < cs.size(); j++) {
				s.preS.get(j).pop += cs.size() - j;
			}
		}

		ArrayList<College> clcs = new ArrayList<College>(cs);
		PLmergeSort(clcs);

		for (int i = 0; i < clcs.size(); i++) {
			for (int j = 0; j < sts.size(); j++) {
				PL.add(new Contract(clcs.get(i).preC.get(j), clcs.get(i)));
			}
		}
	}

	// 全学生の全選好反映、もちろん戦略的操作可能
	void preStuPLSet(ArrayList<Contract> PL, ArrayList<Student> sA, ArrayList<College> cs, ArrayList<Region> rs,
			ArrayList<Student> sts) {
		sA=new ArrayList<Student>();
//		int p=(int)Math.random()*512;
//		sA.add(sts.get(p));
		for(Student s:sts){
			if(Math.random()<siki)
				sA.add(s);
		}
		System.out.println("size: "+sA.size());
		for (int i = 0; i < sA.size(); i++) {
			Student s = sA.get(i);
			for (int j = 0; j < cs.size(); j++) {
				s.preS.get(j).pop += cs.size() - j;
			}
		}

		ArrayList<College> clcs = new ArrayList<College>(cs);
		PLmergeSort(clcs);

		for (int i = 0; i < clcs.size(); i++) {
			for (int j = 0; j < sts.size(); j++) {
				PL.add(new Contract(clcs.get(i).preC.get(j), clcs.get(i)));
			}
		}
	}

	// 謝罪ものの付属品
	void greatSPLSet(ArrayList<Contract> PL, ArrayList<Student> sts, ArrayList<College> cs) {
		ArrayList<College> clcs = new ArrayList<College>(cs);
		for (int i = 0; i < clcs.size(); i++) {
			College c = clcs.get(i);
			if (c.pop > 0) {
				for (int j = 0; j < sts.size(); j++) {
					PL.add(new Contract(c.preC.get(j), c));
				}
			} else
				break;
		}
		for (int i = 0; i < sts.size(); i++) {
			for (int j = 0; j < clcs.size(); j++) {
				College c = clcs.get(j);
				if (c.pop > 0) {
					continue;
				} else
					PL.add(new Contract(c.preC.get(i), c));
			}
		}
	}

	// 学生全員がマッチしたかどうか
	boolean isAllStuMatch(ArrayList<Student> sts) {
		for (int i = 0; i < sts.size(); i++) {
			if (Objects.equals(null, sts.get(i).c))
				return false;
		}
		return true;
	}

	// 学校cが属している地域に下要素的下限があるかどうか
	boolean Ear(College c, ArrayList<Region> rs) {
		for (int i = 0; i < rs.size(); i++) {
			if (rs.get(i).hasC.contains(c))
				if (rs.get(i).ar > 0)
					return true;
		}
		return false;
	}

	// 木のアップデート(詳しくは横尾教授のSDRQ、MSDARQの論文で)
	void update(ArrayList<College> Xc, ArrayList<Integer> Xnum, ArrayList<Region> rs) {
		for (int i = 0; i < Xc.size(); i++) {
			for (int j = 0; j < rs.size(); j++) {
				if (rs.get(j).hasC.contains(Xc.get(i))) {
					rs.get(j).pr = Math.max(0, rs.get(j).pr - Xnum.get(i));
				}
			}
		}
		for (int i = 0; i < rs.size(); i++) {
		}
		for (int i = 0; i < Xc.size(); i++) {
			for (int j = 0; j < rs.size(); j++) {
				if (rs.get(j).hasC.contains(Xc.get(i))) {
					if (rs.get(j).hasC.size() == 1)
						rs.get(j).qr = rs.get(j).qr - Xnum.get(i);
				}
			}
		}
		for (int i = 0; i < Xc.size(); i++) {
			for (int j = 0; j < rs.size(); j++) {
				if (rs.get(j).hasC.contains(Xc.get(i))) {
					if (rs.get(j).hasC.size() > 1)
						rs.get(j).qr = childrenQ(rs.get(j), rs);
				}
			}
		}
		rs.get(0).pr = getRepr(rs.get(0));
		for (int i = 0; i < rs.size(); i++) {
			if (rs.get(i).hasC.size() == 1) {
				rs.get(i).ar = rs.get(i).pr;
				continue;
			}
			rs.get(i).ar = Math.max(0, rs.get(i).pr - (rs.get(i).childrenN[0].pr + rs.get(i).childrenN[1].pr));
		}
	}

	// 地域の学生数計算
	int childrenQ(Region r, ArrayList<Region> rs) {
		int sum = 0;
		for (int i = 0; i < r.hasC.size(); i++) {
			College c = r.hasC.get(i);
			for (int j = 0; j < rs.size(); j++) {
				if (rs.get(j).hasC.contains(c)) {
					if (rs.get(j).hasC.size() == 1)
						sum += rs.get(j).qr;
				}
			}
		}
		return sum;
	}

	// 子ノードの下限の合計より自分の下限が小さい時、自分の下限を修正する
	int getRepr(Region node) {
		if (node.hasC.size() == 1)
			return node.pr;
		if (node.pr < getRepr(node.childrenN[0]) + getRepr(node.childrenN[1]))
			node.pr = getRepr(node.childrenN[0]) + getRepr(node.childrenN[1]);
		return node.pr;
	}

	// 木における深さ計算
	int depth(ArrayList<Region> rs, int i) {
		Region r0 = rs.get(0);
		Region r = rs.get(i);
		double dr = Math.log(r0.hasC.size()) / Math.log(2) - Math.log(r.hasC.size()) / Math.log(2);
		return (int) dr;
	}

	int pSetTree(Region r) {
		if (r.hasC.size() == 1) {
			r.pr = 0;
			return 0;
		}
		int sum = pSetTree(r.childrenN[0]) + pSetTree(r.childrenN[1]) + 4;
		r.pr = sum;
		return sum;
	}

	// 期待最小カウント計算
	int minimumCount(Region r) {
		if (r.hasC.size() == 1) {
			r.er = Math.max(allCollegeStudent(r), r.pr);
			return r.er;
		}
		r.er = Math.max(minimumCount(r.childrenN[0]) + minimumCount(r.childrenN[1]), r.pr);
		return r.er;
	}

	// 妥当な不満を持つ学生の数を計算
	int envy(ArrayList<Student> Sts, int[] cS) {
		int sum = 0;
		for (int i = 0; i < Sts.size(); i++) {
			if (Sts.get(i).maxPre == 0)
				continue;
			int cSum = 0;
			for (int j = Sts.get(i).maxPre - 1; j >= 0; j--) {
				College envyC = Sts.get(i).preS.get(j);
				for (int k = 0; k < envyC.hasS.size(); k++) {
					if (envyC.preC.indexOf(Sts.get(i)) < envyC.preC.indexOf(envyC.hasS.get(k))) {
						//System.out.println(Sts.get(i).number+"->"+envyC.hasS.get(k).number);
						sum++;
						cSum++;
						cS[i]++;
						// System.out.println("number:"+i);
						break;
					}
				}
				if (cSum > 0)
					break;
			}
		}
		// System.out.println("envy:" + sum);
		return sum;
	}
	// 妥当な不満を持つ学生の数を計算
	int envy2(ArrayList<Student> Sts) {
		int sum = 0;
		for (int i = 0; i < Sts.size(); i++) {
			if (Sts.get(i).maxPre == 0)
				continue;
			int cSum = 0;
			for (int j = Sts.get(i).maxPre - 1; j >= 0; j--) {
				College envyC = Sts.get(i).preS.get(j);
				for (int k = 0; k < envyC.hasS.size(); k++) {
					if (envyC.preC.indexOf(Sts.get(i)) < envyC.preC.indexOf(envyC.hasS.get(k))) {
						//System.out.println(Sts.get(i).number+"->"+envyC.hasS.get(k).number);
						sum++;
						cSum++;
//						if(Sts.get(i).group==envyC.hasS.get(k).group)
							envyin++;
	//					else
//							envyout++;
						// System.out.println("number:"+i);
						break;
					}
				}
				if (cSum > 0)
					break;
			}
		}
		// System.out.println("envy:" + sum);
		return sum;
	}

	int speEnvy(ArrayList<Student> Sts, int[] cS) {
		int sum = 0;
		for (int i = 0; i < Sts.size(); i++) {
			if (Sts.get(i).maxPre == 0)
				continue;
			int cSum = 0;
			for (int j = Sts.get(i).maxPre - 1; j >= 0; j--) {
				College envyC = Sts.get(i).preS.get(j);
				for (int k = 0; k < envyC.hasS.size(); k++) {
					if (envyC.preC.indexOf(Sts.get(i)) < envyC.preC.indexOf(envyC.hasS.get(k))) {
						// System.out.println(Sts.get(i).number + "->" +
						// envyC.hasS.get(k).number);
						sum++;
						cSum++;
						cS[i]++;
						break;
					}
				}
				if (cSum > 0)
					break;
			}
		}
		// System.out.println("envy:" + sum);
		return sum;
	}

	int[] allEnvy(ArrayList<Student> Sts, boolean isEnd, int[] cS) {
		int[] chnum = new int[2];
		chnum[0] = 0;
		chnum[1] = 0;
		int max = 0;
		int sum = 0;
		int count = 0;
		int dif = 0;
		int bad = 0;
		int badCount = 0;
		int eqCount = 0;
		int origin = envy(Sts, cS);
		for (int i = 0; i < Sts.size(); i++) {
			Student s1 = Sts.get(i);
			for (int j = i + 1; j < Sts.size(); j++) {
				sum++;
				Student s2 = Sts.get(j);

				College tmpC1 = s1.c;
				College tmpC2 = s2.c;
				int max1 = s1.maxPre;
				int max2 = s2.maxPre;
				ArrayList<Student> tmpLS1 = tmpC1.hasS;
				ArrayList<Student> tmpLS2 = tmpC2.hasS;
				Sts.get(i).c = tmpC2;
				Sts.get(j).c = tmpC1;
				Sts.get(i).maxPre = Sts.get(i).preS.indexOf(Sts.get(i).c);
				Sts.get(j).maxPre = Sts.get(j).preS.indexOf(Sts.get(j).c);
				tmpC1.hasS.remove(s1);
				tmpC2.hasS.remove(s2);
				tmpC1.hasS.add(s2);
				tmpC2.hasS.add(s1);

				int min = envy(Sts, cS) - origin;
				if (min < 0) {
					dif += min;
					count++;
					if (max < -1 * min) {
						max = -1 * min;
						chnum[0] = i;
						chnum[1] = j;
					}
				} else if (min == 0) {
					eqCount++;
				} else {
					bad += min;
					badCount++;
				}
				Sts.get(i).c = tmpC1;
				Sts.get(j).c = tmpC2;
				Sts.get(i).maxPre = max1;
				Sts.get(j).maxPre = max2;
				tmpC1.hasS.remove(s2);
				tmpC2.hasS.remove(s1);
				tmpC1.hasS.add(s1);
				tmpC2.hasS.add(s2);
			}
		}
		System.out.println("maxdif:" + max + ":" + chnum[0] + ";" + chnum[1]);
		if (max == 0)
			isEnd = true;
		return chnum;
	}

	// 空きシートを要求する学生の数を計算
	int waste(ArrayList<Integer> qrs, ArrayList<Integer> prs, ArrayList<Student> Sts, ArrayList<Region> Rs, int[] cS) {
		int sum = 0;
		int[] waste = new int[64];
		for (int i = 0; i < Sts.size(); i++) {
			boolean already = false;
			if (Sts.get(i).maxPre == 0) {
				continue;
			}
			for (int j = Sts.get(i).maxPre - 1; j >= 0 && !already; j--) {
				boolean cNum = false;
				College wasteC = Sts.get(i).preS.get(j);
				wasteC.hasS.add(Sts.get(i));
				Sts.get(i).c.hasS.remove(Sts.get(i));
				for (int k = 0; k < Rs.size(); k++) {
					if (Rs.get(k).hasC.contains(wasteC))
						if (allCollegeStudent(Rs.get(k)) > qrs.get(k)) {
							cNum = true;
						}
				}
				if (!cNum)
					for (int k = 0; k < Rs.size(); k++) {
						if (Rs.get(k).hasC.contains(Sts.get(i).c)) {
							if (allCollegeStudent(Rs.get(k)) < prs.get(k)) {
								cNum = true;
							}
						}
					}
				if (!cNum) {
					already = true;
				}
				wasteC.hasS.remove(Sts.get(i));
				Sts.get(i).c.hasS.add(Sts.get(i));
			}
			if (already) {
				cS[i]++;
				//System.out.println("number:"+i);
				sum++;
			}
		}
		// System.out.println("waste:" + sum);
		return sum;
	}
	// 空きシートを要求する学生の数を計算
	int waste2(ArrayList<Integer> qrs, ArrayList<Integer> prs, ArrayList<Student> Sts, ArrayList<Region> Rs) {
		int sum = 0;
		int[] waste = new int[64];
		for (int i = 0; i < Sts.size(); i++) {
			boolean already = false;
			if (Sts.get(i).maxPre == 0) {
				continue;
			}
			for (int j = Sts.get(i).maxPre - 1; j >= 0 && !already; j--) {
				boolean cNum = false;
				College wasteC = Sts.get(i).preS.get(j);
				wasteC.hasS.add(Sts.get(i));
				Sts.get(i).c.hasS.remove(Sts.get(i));
				for (int k = 0; k < Rs.size(); k++) {
					if (Rs.get(k).hasC.contains(wasteC))
						if (allCollegeStudent(Rs.get(k)) > qrs.get(k)) {
							cNum = true;
						}
				}
				if (!cNum)
					for (int k = 0; k < Rs.size(); k++) {
						if (Rs.get(k).hasC.contains(Sts.get(i).c)) {
							if (allCollegeStudent(Rs.get(k)) < prs.get(k)) {
								cNum = true;
							}
						}
					}
				if (!cNum) {
					already = true;
				}
				wasteC.hasS.remove(Sts.get(i));
				Sts.get(i).c.hasS.add(Sts.get(i));
			}
			if (already) {
				//System.out.println("number:"+i);
				sum++;
			}
		}
		// System.out.println("waste:" + sum);
		return sum;
	}

	// 地域の学生の数を計算
	int allCollegeStudent(Region r) {
		int sum = 0;
		for (int i = 0; i < r.hasC.size(); i++) {
			sum += r.hasC.get(i).hasS.size();
		}
		return sum;
	}

	void showStudentPre(ArrayList<Student> Sts) {
		System.out.println("Student");
		for (int i = 0; i < Sts.size(); i++) {
			System.out.println(i + ":");
			for (int j = 0; j < Student.m; j++) {
				System.out.print(Sts.get(i).preS.get(j).number + ",");
			}
			System.out.println();
		}
	}

	void showCollegePre(ArrayList<College> Cs) {
		System.out.println("");
		System.out.println("College");
		for (int i = 0; i < Cs.size(); i++) {
			System.out.println(i + ":");
			for (int j = 0; j < College.n; j++) {
				System.out.print(Cs.get(i).preC.get(j).number + ",");
			}
			System.out.println();
		}
	}

	double calcSigma(int[] hs, double ave) {
		double sum = 0;
		for (int h : hs) {
			sum += Math.pow(ave - h, 2);
		}
		return Math.sqrt((sum / hs.length));
	}

	void showRegionNum(ArrayList<Region> Rs) {
		for (int i = 0; i < Rs.size(); i++) {
			System.out.println(i + ":" + Rs.get(i).pr);
		}
		for (int i = 0; i < Rs.size(); i++) {
			System.out.println(i + ":" + Rs.get(i).qr);
		}
		for (int i = 0; i < Rs.size(); i++) {
			System.out.println(i + "a:" + Rs.get(i).ar);
		}
	}

	void showDetailMatch(ArrayList<College> Cs) {
		for (int i = 0; i < Cs.size(); i++) {
			System.out.println(Cs.get(i).number);
			for (int j = 0; j < Cs.get(i).hasS.size(); j++) {
				System.out.print(j + ":" + Cs.get(i).hasS.get(j).number + ",");
			}
			System.out.println("");
		}
	}

	void showRoughMatch(ArrayList<College> Cs) {
		String s = "";
		for (int i = 0; i < Cs.size(); i++) {
			s += i + ":" + Cs.get(i).hasS.size() + ", ";
		}
		System.out.println(s);
	}

	// こっから下マージソート(少しでも変わると(昇順、降順など)新しく作ったので数が多く見る必要は特にないです)
	void merge(ArrayList<Student> cs1, ArrayList<Student> cs2, ArrayList<Student> cs, College c) {
		int i = 0, j = 0;
		while (i < cs1.size() || j < cs2.size()) {
			if (j >= cs2.size() || (i < cs1.size() && c.preC.indexOf(cs1.get(i)) < c.preC.indexOf(cs2.get(j)))) {
				cs.set(i + j, cs1.get(i));
				i++;
			} else {
				cs.set(i + j, cs2.get(j));
				j++;
			}
		}
	}

	void mergeSort(int fix, ArrayList<Student> cs, College c) {
		if (cs.size() > 1) {
			ArrayList<Student> fixS = new ArrayList<Student>();
			for (int i = 0; i < fix; i++) {
				fixS.add(cs.get(0));
				cs.remove(0);
			}
			int m = cs.size() / 2;
			int n = cs.size() - m;
			ArrayList<Student> cs1 = new ArrayList<Student>();
			ArrayList<Student> cs2 = new ArrayList<Student>();
			for (int i = 0; i < m; i++)
				cs1.add(cs.get(i));
			for (int i = 0; i < n; i++)
				cs2.add(cs.get(m + i));
			mergeSort(0, cs1, c);
			mergeSort(0, cs2, c);
			merge(cs1, cs2, cs, c);
			for (int i = fix; i > 0; i--) {
				cs.add(0, fixS.get(i - 1));
			}
		}
	}

	void merge(ArrayList<Contract> cs1, ArrayList<Contract> cs2, ArrayList<Contract> cs) {
		int i = 0, j = 0;
		while (i < cs1.size() || j < cs2.size()) {
			if (j >= cs2.size() || (i < cs1.size() && cs1.get(i).getPriority() < cs2.get(j).getPriority())) {
				cs.set(i + j, cs1.get(i));
				i++;
			} else {
				cs.set(i + j, cs2.get(j));
				j++;
			}
		}
	}

	void mergeSort(ArrayList<Contract> cs) {
		if (cs.size() > 1) {
			int m = cs.size() / 2;
			int n = cs.size() - m;
			ArrayList<Contract> cs1 = new ArrayList<Contract>();
			ArrayList<Contract> cs2 = new ArrayList<Contract>();
			for (int i = 0; i < m; i++)
				cs1.add(cs.get(i));
			for (int i = 0; i < n; i++)
				cs2.add(cs.get(m + i));
			mergeSort(cs1);
			mergeSort(cs2);
			merge(cs1, cs2, cs);
		}
	}

	void PLmerge(ArrayList<College> cs1, ArrayList<College> cs2, ArrayList<College> cs) {
		int i = 0, j = 0;
		while (i < cs1.size() || j < cs2.size()) {
			if (j >= cs2.size() || (i < cs1.size() && cs1.get(i).pop > cs2.get(j).pop)
					|| (i < cs1.size() && cs1.get(i).pop == cs2.get(j).pop && cs1.get(i).number < cs2.get(j).number)) {
				cs.set(i + j, cs1.get(i));
				i++;
			} else {
				cs.set(i + j, cs2.get(j));
				j++;
			}
		}
	}

	void PLmergeSort(ArrayList<College> cs) {
		if (cs.size() > 1) {
			int m = cs.size() / 2;
			int n = cs.size() - m;
			ArrayList<College> cs1 = new ArrayList<College>();
			ArrayList<College> cs2 = new ArrayList<College>();
			for (int i = 0; i < m; i++)
				cs1.add(cs.get(i));
			for (int i = 0; i < n; i++)
				cs2.add(cs.get(m + i));
			PLmergeSort(cs1);
			PLmergeSort(cs2);
			PLmerge(cs1, cs2, cs);
		}
	}

	void popmerge(ArrayList<Student> cs1, ArrayList<Student> cs2, ArrayList<Student> cs) {
		int i = 0, j = 0;
		while (i < cs1.size() || j < cs2.size()) {
			if (j >= cs2.size() || (i < cs1.size() && cs1.get(i).pop < cs2.get(j).pop)
					|| (i < cs1.size() && cs1.get(i).pop == cs2.get(j).pop && cs1.get(i).number < cs2.get(j).number)) {
				cs.set(i + j, cs1.get(i));
				i++;
			} else {
				cs.set(i + j, cs2.get(j));
				j++;
			}
		}
	}

	void popmergeSort(ArrayList<Student> cs) {
		if (cs.size() > 1) {
			int m = cs.size() / 2;
			int n = cs.size() - m;
			ArrayList<Student> cs1 = new ArrayList<Student>();
			ArrayList<Student> cs2 = new ArrayList<Student>();
			for (int i = 0; i < m; i++)
				cs1.add(cs.get(i));
			for (int i = 0; i < n; i++)
				cs2.add(cs.get(m + i));
			popmergeSort(cs1);
			popmergeSort(cs2);
			popmerge(cs1, cs2, cs);
		}
	}

	// 謝罪ものの付属品なので気にしないでください。
	void imDA(ArrayList<Student> sts, ArrayList<College> cs, ArrayList<Region> rs, int andNum) {
		ArrayList<Integer> fixNum = new ArrayList<Integer>(); // MSDAQRで使う(固定した前までのマッチングを保存する)
		for (int i = 0; i < cs.size(); i++) { // 各学校の保存生徒数を格納
			fixNum.add(cs.get(i).hasS.size());
		}
		ArrayList<College> Xc = new ArrayList<College>();
		ArrayList<Integer> Xnum = new ArrayList<Integer>();
		while (!(isAllStuMatch(sts))) {// DA
			for (int i = 0; i < cs.size(); i++) {
				int max = cs.get(i).hasS.size();
				for (int j = 0; j < max - fixNum.get(i); j++) {
					cs.get(i).hasS.remove(cs.get(i).hasS.size() - 1);
				}
			}
			Xc = new ArrayList<College>();
			Xnum = new ArrayList<Integer>();

			for (int i = 0; i < sts.size(); i++) {
				Student s = sts.get(i);
				while (true) {
					if (s.preS.get(s.maxPre).forbidden)
						s.maxPre++;
					else
						break;
				}
				s.preS.get(s.maxPre).hasS.add(s);
				s.c = s.preS.get(s.maxPre);
				if (!Xc.contains(s.c)) {
					Xc.add(s.c);
					Xnum.add(1);
				} else {
					Xnum.set(Xc.indexOf(s.c), Xnum.get(Xc.indexOf(s.c)) + 1);
				}
			}
			for (int i = 0; i < cs.size(); i++) {
				College c = cs.get(i);
				int fix = fixNum.get(i);
				mergeSort(fix, c.hasS, c);
				if (c.hasS.size() > andNum) {
					int rem = c.hasS.size() - andNum;
					for (int j = 0; j < rem; j++) {
						c.hasS.get(andNum).c = null;
						c.hasS.get(andNum).maxPre++;
						c.hasS.remove(andNum);
					}
				}
			}
		}
	}
}
