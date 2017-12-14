package b04;
import java.util.ArrayList;

public class Student {
	static int m;// 学校数
	static int sum = 0;// 生徒数(生徒番号に使う)
	int number;// 生徒番号
	ArrayList<College> preS = new ArrayList<College>();// 好みランク
	// ArrayList<Boolean> ispreS = new ArrayList<Boolean>();//
	// 既に応募したか(SDとかでは使わない多分)
	int maxPre = 0;// 現時点での応募できる最高の選好順位
	College c;// マッチする学校
	double u;// 学校からの評価値
	double uS;// 学校からの評価値
	int permit = 0;
	int pop=0;

	int group;

	// コンストラクタ
	public Student() {
		number = sum;
		sum++;
	}
	public Student(double uS) {
		this.uS=uS;
		number = sum;
		sum++;
	}

	// セットアップ(好みまで決める)
	void setUp(double a, ArrayList<College> cs, boolean isImcomplete) {
		permit = 1 + (int) (Math.random() * cs.size());
		ArrayList<College> cst = (ArrayList<College>) cs.clone();
		for (int i = 0; i < m; i++) {
			cst.get(i).u = a * cs.get(i).uC + (1 - a) * Math.random();
			// m次元のベクトルuに評価値と学校のMapを降順で
		}

		mergeSort(cst);
		if (isImcomplete) {
			for (int i = 0; i < permit; i++) {
				preS.add(cst.get(i));
			}
		} else {
			for (int i = 0; i < m; i++) {
				preS.add(cst.get(i));
			}
		}
	}

	void setUp(double a, ArrayList<College> cs) {
		ArrayList<College> cst = (ArrayList<College>) cs.clone();
		for (int i = 0; i < m; i++) {
			cst.get(i).u = a * cs.get(i).uC + (1 - a) * Math.random();
			// m次元のベクトルuに評価値と学校のMapを降順で
		}
		mergeSort(cst);
		for (int i = 0; i < m; i++) {
			preS.add(cst.get(i));
		}
	}

	// こっから下はマージソート
	void merge(ArrayList<College> cs1, ArrayList<College> cs2,
			ArrayList<College> cs) {
		int i = 0, j = 0;
		while (i < cs1.size() || j < cs2.size()) {
			if (j >= cs2.size()
					|| (i < cs1.size() && cs1.get(i).u > cs2.get(j).u)) {
				cs.set(i + j, cs1.get(i));
				i++;
			} else {
				cs.set(i + j, cs2.get(j));
				j++;
			}
		}
	}

	void mergeSort(ArrayList<College> cs) {
		if (cs.size() > 1) {
			int m = cs.size() / 2;
			int n = cs.size() - m;
			ArrayList<College> cs1 = new ArrayList<College>();
			ArrayList<College> cs2 = new ArrayList<College>();
			for (int i = 0; i < m; i++)
				cs1.add(cs.get(i));
			for (int i = 0; i < n; i++)
				cs2.add(cs.get(m + i));
			mergeSort(cs1);
			mergeSort(cs2);
			merge(cs1, cs2, cs);
		}
	}
}
