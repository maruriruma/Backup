package b04;

public class Contract {
	private College c;
	private Student s;
	private int priority;
	private static int sum=0;

	public Contract(Student s,College c){
		this.s=s;
		this.c=c;
		this.priority=sum;
		sum++;
	}
	public Contract(Student s,College c, int p){
		this.s=s;
		this.c=c;
		this.priority=p;
	}
	Student getStudent(){
		return this.s;
	}
	College getCollege(){
		return this.c;
	}
	int getPriority(){
		return this.priority;
	}
	void setPriority(int p){
		this.priority=p;
	}
	boolean contractContain(Contract x){
		if((this.s.equals(x.s))&&(this.c.equals(x.c)))
			return true;
		else
			return false;
	}
}
