#include<iostream>
#include<string>
#include<vector>
#include<math.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using boost::property_tree::ptree;
boost::property_tree::ptree FinalStates;
boost::property_tree::ptree Transitions;
boost::property_tree::ptree StartState;
boost::property_tree::ptree Alphas;
boost::property_tree::ptree StatePointMaps;
boost::property_tree::ptree ControlPointMaps;

using namespace std;
ptree tree;
/////////////////////////////////////////////////////ptree/////////////////////////////////////
void addState(int id, string name);
void addTransition();
void addFinalState();

void initAddFinalState();
void addFinalState(int id, string name);
void finishAddFinalState();
void initAddFinalState();

void initAddTransition();
void addTransition(int idFrom,string nameFrom,int idTo, string nameTo, string input);
void finishAddTransition();

void addStartState(int id, string name);

void initAddAlpha();
void addAlpha(string name);
void finishAddAlpha();

void initStatePointMap();
void addStatePointMap(int id,float x, float y);
void finishAddStatePointMap();

void initAddControlPointMap();
void addControlPointMap(int idFrom,int idTo, float x, float y);
void finishAddControlPointMap();
void init();
void finish();
////////////////////////////////////////////////////////////////////////////////////////
int join[4][4];
string a[4][1000];
typedef struct {                       // Tao cau truc de chua cac phan tu cua tap Closure cua formula
	string phantu[50];
	string phantu_dao[50];
} Closure;
/*
* @brief Tim dau trong bieu thuc string va luu tru cac vi tri dau  vao vecto vtdau
*/
void timDau(string string1, vector<int>&vtdau);
void timDauNgoac(string string1, vector<int>&vtmo, vector<int>&vtdong);// Tim dau ngoac trong bieu thuc string va luu tru cac vi tri dau  vao vecto vtdau
void timClosure(string string1, Closure *Cl, vector<int>&vtdau, int &Clsize, vector<string>&AP); //Tim tap Closure , luu tru trong Cl cua bieu thuc khong co dau ngoac
void timClosurefull(string string1, Closure *Cl, vector<int>&vtdau, vector<int>&vtmo, vector<int>&vtdong, int &Clsize, vector<string>&AP);// Tim closure cua bieu thuc co dau ngoac
void in(string *a, int Clsize);  //In mang a co do lon la Clsize
void in1(vector<string>F, int Clsize){
	for (int i = 0; i < F.size()/Clsize; i++){
		cout<<"\n";
		for (int j = 0; j < Clsize; j++){
			cout << F[i*Clsize + j] << ",";
		}
	}
	cout << '\n';
}
void doinhiphan(int a, int *np, int Clsize);  // Chuyen mot so thap phan a sang so nhi phan co do dai Clsize bit, luu trong mang np
string *doichuoi(int *np, string *str, int Clsize, Closure Cl); //Doi chuoi nhi phan sang cac atom
/*
*@bref : Ham nay kiem tra xem phan tu check co xuat hien trong atom hay khong, neu co thi se tra ve gia tri search khac 0, neu khong thi search =-1
*/
int kiemtraxuathien(string *atom, int atom_size, string check){
	int search = -1;
	for (int i = 0; i < atom_size; i++){
		if (atom[i].compare(check) == 0){
			search = 1;
			break;
		}
	}
	return search;
}
int kiemtraxuathien1(vector<string>atom, string check){
	int search = -1;
	for (int i = 0; i < atom.size(); i++){
		search = atom[i].find(check);
		if (search != -1)
			break;
	}
	return search;
}
int kiemtraxuathien3(vector<string>atom, string check,int start){
	int search = -1;
	for (int i = start; i < atom.size(); i++){
		search = atom[i].find(check);
		if (search != -1)
			break;
	}
	return search;
}
int sosanh2vector(vector<string>atom1, vector < string> atom2){
	int check = 1;
	for (int i = 0; i < atom1.size(); i++){
		if (atom1[i].compare(atom2[i]) != 0)
			check = -1;
	}
	return check;
}

int kiemtraxuathien2(vector<string>atom, vector<string>Q){
	int check = -1;
	vector<string>dem;
	int size1 = Q.size() / atom.size();
	for (int i = 0; i < size1; i++){
		for (int j = 0; j < atom.size(); j++){
			dem.push_back(Q[i*atom.size() + j]);
		}
		if (sosanh2vector(atom, dem) == 1){
			check = 1;
			break;
		}
		dem.clear();
	}
	return check;
}


/*
*@bref : Ham nay them 1 atom vao vecto
*/
void thematomvaovector(string *atom, int Clsize, vector<string>&vector1){
	for (int i = 0; i < Clsize; i++){
		vector1.push_back(atom[i]);
	}
}
void thematomvaovector1(vector<string>atom, vector<string>&vector1){
	for (int i = 0; i < atom.size(); i++){
		vector1.push_back(atom[i]);
	}
}

/*
*@bref : Ham nay lay 1 atom tu tap W va tra lai atom do
*/
vector<string> layatom(int vtatom, vector<string>&vector1, int Clsize){
	vector<string>dem_atom;
	for (int j = 0; j < Clsize; j++){
		dem_atom.push_back(vector1[vtatom*Clsize]);
		vector1.erase(vector1.begin() + vtatom*Clsize);
	}
	return dem_atom;
}
vector<string> timgiaoAP(vector<string>atom, vector<string>AP, int Clsize){
	string giao;
	vector<string>giao1;
	vector<string> dem;
	for (int i = 0; i < AP.size(); i++){
		dem.push_back(AP[i]);
		if (kiemtraxuathien2(dem, atom) == 1){
			giao = giao + AP[i];

		}
		dem.clear();
	}
	
	if (giao.empty() == 1)
		giao1.push_back("O");
	else
		giao1.push_back(giao);
	return giao1;
}
/*
*@bref : Tim vi tri cua U trong 1 string, tra lai mang chua vi tri U, phan tu dau tien cua mang chi so luong cua mang
		Bat dau luu vi tri tu phan tu [1]
*/
int *timU(string str){
	int *a = new int;
	int count = 0;
	for (int i = 0; i < str.size(); i++){
		a[0] = count;
		if (str[i] == 'U'){
			count++;
			a[count] = i;
		}
	}
	return a;
}
int *timngoac(string str){
	int *a = new int;
	int count = 0;
	for (int i = 0; i < str.size(); i++){
		a[0] = count;
		if (str[i] == 'U'){
			count++;
			a[count] = i;
		}
	}
	return a;
}

void check_line_6vs7(vector<string> &F_p1Up2,Closure Cl,int Clsize,vector<string> atom_lay);


void tinhJoin(int m, int n){		// tham so la mang anpha i, m, n
	for (int i = 1; i <= m; i++){
		for (int j = 1; j <= n; j++){
			if (a[i][j][0] == 'X') {	//Tim thay phep Xfi
				string fi;
				fi = a[i][j].substr(1, a[i][j].length() - 1);			// Neu Xfi thi lay fi
				if (fi[0] == '(') fi = fi.substr(1, fi.length() - 2);  // Neu X(fi) thi xoa dau nguoac de lay fi
				for (int x = 1; x <= m; x++) {
					bool ok = false;
					for (int y = 1; y <= n; y++)
					if (fi == a[x][y]) {			// Duyet phan tu cua anpha[x] tim fi. neu tim thay thi ok=true, break
						ok = true;
						break;
					}
					if (ok == false) join[i][x] = 0;	// Neu fi ko thuoc anpha[x] thi ko noi i voi x
				}
			}
		}
	}

	string a1, a2, tmp, a1Ua2;
	int chicky1 = 0;//bien the hien cho (aUb)Uc va ta dang can xet den U thu 2
	int chicky2 = 0;// biens the hien cho vi chi cua ) ben tren neu co vaf neu ko co la -1
	tmp = a[1][n];
	if (tmp[0] == '!') tmp = tmp.substr(1, tmp.length() - 1);
	if (tmp[0] == '(') {
		chicky1 = 1;
		chicky2 = tmp.find(')');
	}
	int indexU = tmp.find('U', chicky2);
	int chicky3 = 0;
	if (tmp[indexU + 1] == '(') chicky3 = 1;
	if (indexU != -1){
		if (chicky1 == 1) a1 = tmp.substr(1, chicky2 - 1);
		else a1 = tmp.substr(0, indexU);
		a2 = tmp.substr(indexU + 1 + chicky3, tmp.length() - indexU - 3 * chicky3);
		a1Ua2 = a1 + 'U' + a2;
		bool check = false;
		bool ok1 = false;
		bool ok2 = false;
		if (a1Ua2 == a[1][n]) check = true;
		for (int j = 1; j <= n; j++){
			if (a[1][j] == a1) ok1 = true;
			if (a[1][j] == a2) ok2 = true;
		}

		if (check == true){
			if (ok2 == true) {
				join[1][2] = join[1][2] && 1;// neu aUb va b deu thuoc atom _lay thi join=1
			}
			else {
				if (ok1 == true) { // neu khong 
					int chicky = 0;
					for (int x = 1; x <= n; x++)
					if ((a[2][x] == a1Ua2)) chicky = 1;
					if (chicky == 0) join[1][2] = join[1][2] && 0;
				}
				else join[1][2] = join[1][2] && 0;
			}
		}
		else {// check ==false
			if (ok2 == true) {
				join[1][2] = join[1][2] && 0;
			}
			if (ok1 == true){
				for (int j = 1; j <= n; j++){
					if (a[2][j] == a1Ua2) join[1][2] = join[1][2] && 0;
				}
			}
		}
	}
}
int check1112(int Clsize, vector<string> atom_lay, vector<string> beta){
	int m ,n;
	m=3;
	n = Clsize;
	for (int count=1; count<= Clsize; count++){
		a[1][count] = atom_lay[count-1];
		a[2][count] = beta[count-1];
	}
	for (int i=1;i<=m;i++)
		for (int j=1;j<=m;j++) join[i][j] = 1; // Khoi tao mang Join
			tinhJoin(m,n);	// goi Ham
	if(join[1][2] == 1) return 1;
	else return 0;
}
int timchiso(string a,string *b,int size_pat){
	int k;
	for(int i=0;i<size_pat;i++){
		if (a.compare(b[i])==0)
			k=i;
	}
	return k;
}

int main()
{
	string formula;
	vector<int>vtmo;					// Vi tri dau mo ngoac
	vector<int>vtdong;					// Vi tri dau dong ngoac
	vector<int>vtdau;					// Chuoi chua cac vi tri co dau
	vector<string>AP;					// Bang chu cai AP
	Closure Cl;
	int Clsize = 0;

	// Nhap bieu thuc formula
	cout << "\nNhap vao bieu thuc:" << endl;
	cout << "Chu y: Cac dau duoc coi la co nghia neu chung la !(phu dinh), U(until), X(next), A(and)" << endl;
	getline(cin, formula, '\n');
	timDauNgoac(formula, vtmo, vtdong);
	while (vtmo.size() != vtdong.size()){
		cout << "Thieu dau ngoac";
		cout << "\nNhap vao bieu thuc:" << endl;
		cout << "Chu y: Cac dau duoc coi la co nghia neu chung la !(phu dinh), U(until), X(next), A(and)" << endl;
		getline(cin, formula, '\n');
	}

	//Tim Dau trong bieu thuc roi di tim tap Closure
	timDau(formula, vtdau);
	if (vtmo.size() == 0)
	{
		timClosure(formula, &Cl, vtdau, Clsize, AP);
	}
	else
		timClosurefull(formula, &Cl, vtdau, vtmo, vtdong, Clsize, AP);

	//Modify Cl.phantu_dao
	for (int i = 0; i < Clsize; i++){
		if (Cl.phantu[i][0] == '!'){
			for (int j = i+1; j < Clsize; j++){
				Cl.phantu_dao[j - 1] = Cl.phantu_dao[j];
				Cl.phantu[j - 1] = Cl.phantu[j];
			}
			Clsize--;
		}
	}
	//In tap closure
	cout << "Tap Closure: ";
	for (int i = 0; i < Clsize; i++){
		cout << Cl.phantu[i] << ",";
		cout << Cl.phantu_dao[i] << ",";
	}



	//Tim atom phu - chua kiem tra dieu kien;
	int size = pow(2, Clsize);
	int *np = new int[Clsize];
	string *str = new string[size*Clsize];
	string **pstr = new string *[size];
	for (int i = 0; i < size; i++){
		doinhiphan(i, np, Clsize);   //tra ve vi tri mang nhi phan
		pstr[i] = doichuoi(np, str + i*Clsize, Clsize, Cl);
	}
	
	// Kiem  tra dieu kien de tim atom tu tap pstr


	//Bat dau khai trien thuat toan
	// Nen kiem tra dieu kien 12 truoc, neu at nao khong thoa man thi bo luon

	string *at = new string[size*Clsize];
	string **pat = new string *[size];
	int size_pat = 0;


	//Tam thoi
	at = str;
	pat = pstr;
	size_pat = size;
	string at_tree[50];


	cout << "\nCo " << size_pat << " atom\n";
	for (int k = 0; k < size_pat; k++){
		for (int i = 0; i < Clsize; i++){
			cout << pat[k][i] << ",";
			at_tree[k]+=pat[k][i];
			at_tree[k]+=",";

		}
		cout << '\n';
	}

	// Buoc 1 : Gan Qo<- at; Q<-0; sigma <- 0;
	vector<string>Qo;
	vector<string>Q;
	vector<string>sigma;
	for (int i = 0; i < size_pat; i++){
		if (kiemtraxuathien(pat[i], Clsize, formula) == 1)
			thematomvaovector(pat[i], Clsize, Qo);
	}
	//in1(Qo,Clsize);  //--> Checked
	// Buoc 2 : Tao W -<Qo
	vector<string>W(Qo);
	//in1(W, Clsize);  //-->Checked
	vector<string>F_p1Up2;
	//Buoc 3: Kiem tr W!=0;
	//
	int bo = 0;
	//while (bo < 1){ -->de check
	while (!(W.empty())){
		//Buoc 4 : Pick atom alpha khoi W
		vector<string>atom_lay;
		atom_lay.clear();
		atom_lay = layatom(0, W, Clsize);
		vector<string>dem_giaoAP = timgiaoAP(atom_lay, AP, Clsize);
		//in1(atom_lay, Clsize);  --> Checked
		//in1(W, Clsize);  --> Checked
		//Buoc 5: Them atom vua lay vao Q
		if (kiemtraxuathien2(atom_lay, Q) == -1)
			thematomvaovector1(atom_lay, Q);
		//in1(Q, Clsize); cout << "t";--> Checked
		
		//Buoc6 + 7;
		check_line_6vs7(F_p1Up2,Cl,Clsize,atom_lay);
        //cout<<"\nF cua phep U la:";

		//Buoc 8 quet tap at
		vector<string>beta;
		//string *beta = new string[Clsize];
		for (int i = 0; i < size_pat; i++){
			for (int j = 0; j < Clsize; j++){
				beta.push_back(pat[i][j]);
			}
			// Buoc 9
			//Buoc 10
			//cout<<"\natom_lay la"<<endl;
			//in1(atom_lay,Clsize);
			//cout<<"\nbeta"<<endl;
			//in1(beta,Clsize);
			//cout<<"\ncheck1112la: "<<check1112(Clsize, atom_lay,beta)<<"\n"<<endl;
		//	std::cin.get();
			if (check1112(Clsize, atom_lay,beta)==1){
				vector<string>dem_chuyen;

				thematomvaovector1(atom_lay, dem_chuyen);
				thematomvaovector1(beta, dem_chuyen);
				thematomvaovector1(dem_giaoAP, dem_chuyen);
				if (kiemtraxuathien2(dem_chuyen, sigma) == -1)
				{
					thematomvaovector1(dem_chuyen, sigma);
				}
				dem_chuyen.clear();

				//Buoc 11
				if (kiemtraxuathien2(beta, Q) == -1){
					thematomvaovector1(beta, W);
				}
				beta.clear();
			}
			//beta.clear();
			
		}
		//in1(sigma, 2 * Clsize + 1); //-->Checked
		//in1(W, Clsize);//--> Checked
		//bo++;--> de check
		dem_giaoAP.clear();
	}
	cout << "\n Tap cac trang thai Q:"; in1(Q, Clsize);
	cout << "\n Tap cac trang thai bat dau Qo:"; in1(Qo, Clsize);
	cout << "\n Tap cac trang thai ket thuc F:"; in1(F_p1Up2, Clsize);
	cout << "\n Bang chu cai AP:"; in1(AP, AP.size());
	cout << "\n Tap cac ham chuyen sigma:"; in1(sigma, 2 * Clsize + 1);


cin.get();
////// for tree
int Qo_count = Qo.size()/Clsize;
string Qo_tree[50];
for(int i=0;i<Qo_count;i++){
	for(int j=0;j<Clsize;j++){
		Qo_tree[i]+=Qo[i*Clsize+j];
		Qo_tree[i]+=",";

	}
}

int F_count = F_p1Up2.size()/Clsize;
string F_tree[50];
for(int i=0;i<F_count;i++){
	for(int j=0;j<Clsize;j++){
		F_tree[i]+=F_p1Up2[i*Clsize+j];
		F_tree[i]+=",";
	}
}

int sigma_count = sigma.size()/(Clsize*2+1);


	

	///////////////////////////// buil ptree and convert to the jflap/////////////////
	init();
	// add state
	for (int i = 0; i<size_pat; i++){
		addState(i,at_tree[i]);
	}
	addState(size_pat,"init");
	// add final state
	initAddFinalState();
    for (int i=0;i<F_count;i++){
    	addFinalState(timchiso(F_tree[i],at_tree,size_pat),F_tree[i]);
    }
    finishAddFinalState();
    //add transition
    //addStartState(timchiso(Qo_tree[0],at_tree,size_pat),Qo_tree[0] );

    // add transition

    initAddTransition();
    for(int i=0;i<sigma_count;i++){
    	string qo,q1;
    	for(int j=0;j<Clsize;j++){
    		qo+=sigma[i*(2*Clsize+1)+j];
    		qo+=",";
    		q1+=sigma[i*(2*Clsize+1)+Clsize+j];
    		q1+=",";
    	}
    	addTransition(timchiso(qo,at_tree,size_pat),qo,timchiso(q1,at_tree,size_pat),q1,sigma[(i+1)*(2*Clsize+1)-1]);
    }
 
    for(int k=0;k<Qo_count;k++){
    		addTransition(size_pat,"init",timchiso(Qo_tree[k],at_tree,size_pat),Qo_tree[k],"lamda");
    	
    }
    //addTransition(0,"q0",1,"q1","a");
    finishAddTransition();
    // add start state
    //addStartState(0,"q0");
    addStartState(size_pat,"init");
    
    // add the alpha
     initAddAlpha();
     for(int i=0;i<AP.size();i++){
     	addAlpha(AP[i]);
     }
    //addAlpha("a");
    finishAddAlpha();
    // add state point map

    //ma tran vi tri
    int vitri[2][50];

    initStatePointMap();
    for (int i=0 ;i<size_pat;i++){
    	int x = 50+5*i;
    	int y=50+i*i;
    	vitri[0][i]=x;
    	vitri[1][i]=y;
    	addStatePointMap(i,x,y);
    }
    addStatePointMap(size_pat,25,25);
    finishAddStatePointMap();
    //
    initAddControlPointMap();
    for(int i=0;i<sigma_count;i++){
    	string qo,q1;
    	for(int j=0;j<Clsize;j++){
    		qo+=sigma[i*(2*Clsize+1)+j];
    		qo+=",";
    		q1+=sigma[i*(2*Clsize+1)+Clsize+j];
    		q1+=",";
    	}
    	if (timchiso(q1,at_tree,size_pat)==timchiso(qo,at_tree,size_pat)){
    		addControlPointMap(timchiso(qo,at_tree,size_pat),timchiso(q1,at_tree,size_pat),0.5*(vitri[0][timchiso(q1,at_tree,size_pat)]+vitri[0][timchiso(q1,at_tree,size_pat)])+60,0.5*(vitri[1][timchiso(q1,at_tree,size_pat)]+vitri[1][timchiso(q1,at_tree,size_pat)]));
    	}
    	else
    		addControlPointMap(timchiso(qo,at_tree,size_pat),timchiso(q1,at_tree,size_pat),0.5*(vitri[0][timchiso(q1,at_tree,size_pat)]+vitri[0][timchiso(q1,at_tree,size_pat)]),0.5*(vitri[1][timchiso(q1,at_tree,size_pat)]+vitri[1][timchiso(q1,at_tree,size_pat)]));

    	
    }
    
     for(int i=0;i<Qo_count;i++){
    	addControlPointMap(size_pat,timchiso(Qo_tree[i],at_tree,size_pat),0.5*(25+vitri[0][timchiso(Qo_tree[i],at_tree,size_pat)]),0.5*(25+vitri[1][timchiso(Qo_tree[i],at_tree,size_pat)]));
    
    }
   
    //addControlPointMap(0,1,50.0,150.0);
    finishAddControlPointMap();
    //finish 
    finish();
	//////////////////////////////////////////////////////////////////////////////////
	return 0;

}
void timDau(string string1, vector<int>&vtdau){
	for (int i = 0; i <= string1.size(); i++) {  // So sanh - tim cac dau
		switch (string1[i]){
		case 'U':
			vtdau.push_back(i);      // Them dau vao chuoi chua dau
			break;
		case 'A':
			vtdau.push_back(i);      // Them dau vao chuoi chua dau
			break;
		case '!':
			vtdau.push_back(i);      // Them dau vao chuoi chua dau
			break;
		case 'X':
			vtdau.push_back(i);      // Them dau vao chuoi chua dau
			break;
		}
	}

}
void timDauNgoac(string string1, vector<int>&vtmo, vector<int>&vtdong){
	for (int i = 0; i <= string1.size(); i++) {  // So sanh - tim cac dau
		switch (string1[i]){
		case '(':
			vtmo.push_back(i);      // Them dau vao chuoi chua dau
			break;
		case ')':
			vtdong.push_back(i);      // Them dau vao chuoi chua dau
			break;
		}
	}
}
void timClosure(string string1, Closure *Cl, vector<int>&vtdau, int &Clsize, vector<string>&AP){
	string closure;
	string dem;
	string phu_dem;

	//Mo rong mang dau
	vector<int>vtdau2 = { -1 };
	for (int i = 1; i <= vtdau.size(); i++){
		vtdau2.push_back(vtdau[i - 1]);
	}
	vtdau2.push_back(string1.size());

	//Bat dau duyet va cat bieu thuc
	for (int i = 0; i <= vtdau.size(); i++){
		int j = i + 1;
		dem = string1.substr(vtdau2[i] + 1, vtdau2[j] - vtdau2[i] - 1);  //Lay ra phan tu
		int check = 1, t;   //Kiem tra xem phan tu do da duoc ke chua
		for (t = 0; t <= Clsize; t++){
			if ((dem.compare((*Cl).phantu[t]) == 0) || (dem.compare(('(' + (*Cl).phantu[t] + ')')) == 0)){
				check = 0;
			}
		}
		if (check != 0){
			(*Cl).phantu[Clsize] = dem; //cout << '\n' << (*Cl).phantu[Clsize];
			(*Cl).phantu_dao[Clsize] = '!' + dem; //cout << '\n' << (*Cl).phantu_dao[Clsize];
			Clsize++;
			AP.push_back(dem);
		}
		dem.clear(); phu_dem.clear();
	}
	int i = 0;
	for (int j = i + 1; j <= vtdau.size() + 1; j++){
		dem = string1.substr(vtdau2[i] + 1, vtdau2[j] - vtdau2[i] - 1);

		//Lay ra phan tu
		int check = 1, t;   //Kiem tra xem phan tu do da duoc ke chua
		for (t = 0; t <= Clsize; t++){
			if ((dem.compare((*Cl).phantu[t]) == 0) || (dem.compare(('(' + (*Cl).phantu[t] + ')')) == 0)){
				check = 0;
			}
		}
		if (check != 0){
			(*Cl).phantu[Clsize] = dem; //cout << '\n' << (*Cl).phantu[Clsize];
			(*Cl).phantu_dao[Clsize] = '!' + dem;// cout << '\n' << (*Cl).phantu_dao[Clsize];
			Clsize++;
		}
		dem.clear(); phu_dem.clear();
	}
}
void timClosurefull(string string1, Closure *Cl, vector<int>&vtdau, vector<int>&vtmo, vector<int>&vtdong, int &Clsize, vector<string>&AP){
	//vtmo1,vtdong1,sl_mo,sl_dong
	string dem;
	vector<int>dem_vtdau;
	//cout << sl_mo - 1;
	while (vtmo.size() > 0){
		for (int i = 0; i <= vtmo.size() - 1; i++){
			int find;
			for (int j = 0; j <= vtdong.size() - 1; j++){
				find = 0;
				if ((vtmo.size() == 1) || ((vtmo[i + 1] > vtdong[j]) && (vtmo[i] < vtdong[j]))){ //Tim cap dong mo ngoac
					dem = string1.substr(vtmo[i] + 1, vtdong[j] - vtmo[i] - 1); //Lay doan trong ngoac
					vector<int>vtxoa;
					for (int k = 0; k < vtdau.size(); k++){
						if ((vtdau[k]> vtmo[i]) && (vtdau[k] < vtdong[j])){
							vtxoa.push_back(k);
						}
					}
					for (int i = 0; i < vtxoa.size(); i++)
						vtdau.erase(vtdau.begin() + vtxoa[0]);
					timDau(dem, dem_vtdau);
					timClosure(dem, Cl, dem_vtdau, Clsize, AP);
					dem.clear();
					dem_vtdau.clear();
					vtmo.erase(vtmo.begin() + i);
					vtdong.erase(vtdong.begin() + j);
					find = 1;
					break;
				}
			}
			if (find = 1)
				break;
		}
	}
	timClosure(string1, Cl, vtdau, Clsize, AP);
}
void in(string *a, int Clsize){
	for (int i = 0; i < Clsize; i++){
		cout << a[i] << ',';
	}
	cout << '\n';
}
void doinhiphan(int a, int *np, int Clsize){
	for (int i = 0; i<Clsize; i++){
		np[i] = 0;
	}
	for (int i = Clsize - 1; a > 0; i--){
		np[i] = a % 2;
		a = a / 2;
	}
}
string *doichuoi(int *np, string *str, int Clsize, Closure Cl){
	for (int i = 0; i <Clsize; i++){
		if (np[i] == 1){
			str[i] = Cl.phantu[i];
		}
		else
			str[i] = Cl.phantu_dao[i];
	}
	return str;
}
void check_line_6vs7(vector<string> &F_p1Up2,Closure Cl,int Clsize,vector<string> atom_lay){
    int count1,count2;
    for(count1 =0;count1< Clsize;count1++){
        vector<string> a;
        int b;
        int c;
        b= Cl.phantu[count1].find('U');
        if(b!=-1){
            if (Cl.phantu[count1].find('U',b+1)==-1)
                    a.push_back(Cl.phantu[count1].substr(Cl.phantu[count1].find('U')+1,Cl.phantu[count1].size()-Cl.phantu[count1].find('U'))) ;
            else a.push_back(Cl.phantu[count1].substr(b+2,Cl.phantu[count1].size()-b-3));
            if ((kiemtraxuathien2(a,atom_lay) !=-1)||( kiemtraxuathien1(atom_lay,Cl.phantu[count1]) == 1)) {
                    if (kiemtraxuathien2(atom_lay,F_p1Up2) ==-1)thematomvaovector1(atom_lay, F_p1Up2);
            }
            a.clear();
            }
        }
}
/////////////////////////////////////////////ptree///////////////////////////////////////////
void init(){
    tree.put("structure.<xmlattr>.type", "editor_panel");
    tree.put("structure.structure.<xmlattr>.type","transition_graph");
    tree.put("structure.structure.structure.<xmlattr>.mode","Default mode");
    tree.put("structure.structure.structure.<xmlattr>.type","fsa");
    tree.put("structure.structure.structure.structure.<xmlattr>.type","state_set");

}
void finish(){
    tree.put("structure.state_label_map","");
    tree.put("structure.note_map","");
    boost::property_tree::xml_writer_settings<std::string> settings('\t',1);
    boost::property_tree::write_xml("./output.jflap",tree,std::locale(), settings);

}
void addState(int id, string name){

    boost::property_tree::ptree ptr;
    ptr.put("name", name);
    ptr.put("id", id);
    tree.add_child("structure.structure.structure.structure.state", ptr);
}

void initAddFinalState(){

    FinalStates.put("<xmlattr>.type","final_states");
}
void addFinalState(int id,string name){
    boost::property_tree::ptree ptr;
    ptr.put("name", name);
    ptr.put("id", id);
    FinalStates.add_child("state", ptr);
}
void finishAddFinalState(){
    tree.add_child("structure.structure.structure.structure",FinalStates);
}

void initAddTransition(){

    Transitions.put("<xmlattr>.type","transition_set");
}
void addTransition(int idFrom,string nameFrom,int idTo, string nameTo, string input){
    boost::property_tree::ptree ptrAdd,ptrFrom,ptrTo;
    ptrAdd.put("input", input);
    ptrFrom.put("name",nameFrom);
    ptrFrom.put("id",idFrom);
    ptrTo.put("name",nameTo);
    ptrTo.put("id",idTo);
    ptrAdd.add_child("from",ptrFrom);
    ptrAdd.add_child("to",ptrTo);
    Transitions.add_child("fsa_trans", ptrAdd);
}
void finishAddTransition(){
    tree.add_child("structure.structure.structure.structure",Transitions);
}

void addStartState(int id, string name){
    StartState.put("<xmlattr>.type","start_state");
    boost::property_tree::ptree ptr;
    ptr.put("name", name);
    ptr.put("id", id);
    StartState.add_child("state", ptr);
    tree.add_child("structure.structure.structure.structure",StartState);
}

void initAddAlpha(){

    Alphas.put("<xmlattr>.type","input_alph");
}
void addAlpha(string name){
    Alphas.add("name", name);
}
void finishAddAlpha(){
    tree.add_child("structure.structure.structure.structure",Alphas);
}

void initStatePointMap(){
}
void addStatePointMap(int id,float x, float y){
    boost::property_tree::ptree ptrPoint,ptrStatePoint;
    ptrStatePoint.put("state",id);

    ptrPoint.put("x",x);
    ptrPoint.put("y",y);
    ptrStatePoint.add_child("point",ptrPoint);
    StatePointMaps.add_child("state_point", ptrStatePoint);
}
void finishAddStatePointMap(){
    tree.add_child("structure.structure.state_point_map",StatePointMaps);
}

void initAddControlPointMap(){

    Transitions.put("<xmlattr>.typr","transition_set");
}
void addControlPointMap(int idFrom,int idTo, float x, float y){
    boost::property_tree::ptree ptrControlPoint, ptrPoint;
    ptrControlPoint.put("from",idFrom);
    ptrControlPoint.put("to",idTo);
    ptrPoint.put("x",x);
    ptrPoint.put("y",y);
    ptrControlPoint.add_child("point",ptrPoint);
    ControlPointMaps.add_child("control_point", ptrControlPoint);
}
void finishAddControlPointMap(){
    tree.add_child("structure.structure.control_point_map",ControlPointMaps);
}
///////////////////////////////////////////////////////////////////////////////////////////////

