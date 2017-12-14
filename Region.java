package b04;
import java.util.ArrayList;

public class Region implements Cloneable {
	Region parentN;// 親ノード
	Region[] childrenN = new Region[2];// 子ノード
	ArrayList<College> hasC = new ArrayList<College>();// 属している学校
	int pr;
	int qr;
	int ar;
	static int aR;
	static int aN;
	static int reN = 4 * 2 - 1;
	static int moreRe;
	static int lessRe;
	int er;

	int nest;

	// コンストラクタ
	public Region(Region pN, ArrayList<College> cs,int q, int n) {
		parentN = pN;
		hasC = (ArrayList<College>) cs.clone();
		ArrayList<College> cs2 = (ArrayList<College>) cs.clone();
		nest=n;
		if (cs2.size() == 1) {
			qr = q;
			pr = 0;
			ar = 0;
			double rand = Math.random();
			if (moreRe == 0) {
				pr = aN - 1;
				ar = aN - 1;
				reN--;
				aR -= ar;
				lessRe--;
			} else if (lessRe == 0) {
				pr = aN;
				ar = aN;
				reN--;
				aR -= ar;
				moreRe--;
			} else if (rand > 0.5) {
				pr = aN;
				ar = aN;
				reN--;
				aR -= ar;
				moreRe--;
			} else {
				pr = aN - 1;
				ar = aN - 1;
				reN--;
				aR -= ar;
				lessRe--;
			}
			Region[] cN = {};
			childrenN = cN;
			return;
		}
		ArrayList<College> cs1 = new ArrayList<College>();
		for (int i = 0; i < (cs.size() / 2); i++) {
			cs1.add(cs2.get(0));
			cs2.remove(0);
		}
		childrenN[0] = new Region(this, cs1,q,n+1);
		childrenN[1] = new Region(this, cs2,q,n+1);
		qr = childrenN[0].qr + childrenN[1].qr;
		pr = childrenN[0].pr + childrenN[1].pr;
		ar = 0;
		double rand = Math.random();
		if (moreRe == 0) {
			pr += aN - 1;
			ar = aN - 1;
			reN--;
			aR -= ar;
			lessRe--;
		} else if (lessRe == 0) {
			pr += aN;
			ar = aN;
			reN--;
			aR -= ar;
			moreRe--;
		} else if (rand > 0.5) {
			pr += aN;
			ar = aN;
			reN--;
			aR -= ar;
			moreRe--;
		} else {
			pr += aN - 1;
			ar = aN - 1;
			reN--;
			aR -= ar;
			lessRe--;
		}
	}
}

