package b04;

import java.util.ArrayList;

public class College {
	static int n;// 生徒数
	static int sum=0;//生徒数(生徒番号に使う)
	int number;//生徒番号
	ArrayList<Student> preC = new ArrayList<Student>();// 好みランク
	int maxPre = 0;// 現時点での応募できる最高の選好順位
	double uC;//世間一般の自校に対する評価値
	double u;//世間と学生個人の評価を合わせたもの
	ArrayList<Student> hasS = new ArrayList<Student>();//マッチングしてる学生
	boolean forbidden=false;
	int pop=0;
	int andNum=Integer.MAX_VALUE;

	// コンストラクタ
	public College(double ut) {
		this.uC=ut;
		number=sum;
		sum++;
	}

	//セットアップ(好みの学生まで)
	void setUp(ArrayList<Student> sts){
		ArrayList<Student> subSts=(ArrayList<Student>) sts.clone();
		for(int i=0;i<n;i++){
			subSts.get(i).u=Math.random();
		}
		mergeSort(subSts);
		for (int i = 0; i < n; i++) {
			preC.add(subSts.get(i));
		}
	}
	void setUp(double b,ArrayList<Student> sts){
		ArrayList<Student> subSts=(ArrayList<Student>) sts.clone();
		for(int i=0;i<n;i++){
			subSts.get(i).u=b*subSts.get(i).uS+(1-b)*Math.random();
		}
		mergeSort(subSts);
		for (int i = 0; i < n; i++) {
			preC.add(subSts.get(i));
		}
	}
	void setUp(ArrayList<Student> sts,int popS){
		pop=popS;
		ArrayList<Student> subSts=(ArrayList<Student>) sts.clone();
		for(int i=0;i<n;i++){
			subSts.get(i).u=Math.random();
			mergeSort(subSts);
		}
		for (int i = 0; i < n; i++) {
			preC.add(subSts.get(i));
		}
	}

	//こっから下はマージソート
	void merge(ArrayList<Student> cs1, ArrayList<Student> cs2, ArrayList<Student> cs) {
		int i = 0, j = 0;
		while (i < cs1.size() || j < cs2.size()) {
			if (j >= cs2.size() || (i < cs1.size() && cs1.get(i).u > cs2.get(j).u)) {
				cs.set(i+j, cs1.get(i));
				i++;
			} else {
				cs.set(i+j, cs2.get(j));
				j++;
			}
		}
	}
	void mergeSort(ArrayList<Student> cs) {
		if (cs.size() > 1) {
			int m = cs.size() / 2;
			int n = cs.size() - m;
			ArrayList<Student> cs1 = new ArrayList<Student>();
			ArrayList<Student> cs2 = new ArrayList<Student>();
			for (int i = 0; i < m; i++)
				cs1.add(cs.get(i));
			for (int i = 0; i < n; i++)
				cs2.add(cs.get(m+i));
			mergeSort(cs1);
			mergeSort(cs2);
			merge(cs1, cs2, cs);
		}
	}
}
