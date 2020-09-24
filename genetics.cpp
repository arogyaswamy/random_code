// Keshav Arogyaswamy
// This program accepts a DNA sequence and outputs the corresponding polypeptide sequence
// If requested, the program will find any Open Reading Frames that may be present in the
// sequence, including all 6 frames (both strands).

#include<iostream>
#include<string>
#include<cctype>
using namespace std;

const string peptide[22] =		{"Ala", "Arg", "Asn", "Asp", "Cys",
								 "Glu", "Gln", "Gly", "His", "Ile",
								 "Leu", "Lys", "Met", "Phe", "Pro",
								 "Ser", "Thr", "Trp", "Tyr", "Val",
								 "STP", "ERR"};

const string one_letter[22]=	{"A",	 "R",	"N",   "D",   "C",
								 "E",	 "Q",	"G",   "H",   "I",
								 "L",	 "K",	"M",   "F",	  "P", 
								 "S",	 "T",	"W",   "Y",   "V",
								 "#",	 "X"};
const int size_of_arrays = 22;
/*
   0: Alanine			5: Glutamic Acid	10: Leucine			15: Serine
   1: Arginine			6: Glutamine		11: Lysine			16: Threonine
   2: Asparigine		7: Glycine			12: Methionine		17: Tryptophan
   3: Aspartic Acid		8: Histidine		13: Phenylalanine	18: Tyrosine
   4: Cysteine			9: Isoleucine		14: Proline			19: Valine

									20: STOP
									21: ERROR

 */

string ORF_checker(string seqn, int start_site, int& end_site, bool& is_min, int min);
string one_letter_conversion(string AA_seqn);
string get_fancy_amino(string codon);
string get_simple_amino(string codon);
string DNA_complement(string template_strand);
void enter_sequence();
void get_output_type();
void ORF_finder(string seqn, int minimum = 0);
void print_polypeptide(string seqn);
void print_pretty_DNA(string seqn);
void print_nt(string seqn, string comp = "NONE");

char output_type;
int length;
int protein_count = 0;
string sequence = "";
string complement = "";

int main(){


	output_type = '0'; //Safeguard to make sure the user has to choose an output type
	const int program_success_code = 0;

	cout << 
		"       _.--\'\"`\'--._    _.--\'\"`\'--._    _.--\'\"`\'--._		\n" << 
		"    \'-:`.\'|`|\"\':-.  \'-:`.\'|`|\"\':-.  \'-:`.\'|`|\"\':-.  '.	\n" <<
		"   .  \'.  | |  | |\'.  \'.  | |  | |\'.  \'.  | |  | |\'.		    \n" <<
		"    \'.  \'.| |  | |  \'.  \'.| |  | |  \'.  \'.| |  | |  \'.		\n" <<
		"      \'._ `.:_ | :_.\' \'._ `.:_ | :_.\' \'._ `.:_ | :_.\'		\n" <<
		"         `-..,..-\'       `-..,..-\'       `-..,..-\'				\n";

	cout << "WELCOME TO KESHAV'S GENETIC ANALYSIS PROGRAM!\n\n";

//Menu
	char menu_option = 'x';
	string dummy;

menu_start:
	cout << "\n\n\n\n\n";
	cout << "Please make a valid selection from the list:\n"
		 << "(Note: You must have an output type selected before printing sequences.)\n"
		 << "0 = Exit		1 = Enter Sequence	2 = Print Sequence/Complement\n"
		 << "3 = Find ORFs		4 = Get AA Seq'n	5 = Select Output Type		 \n";
	cin >> menu_option;


	if(menu_option == '0'){ //Exit
		cout << "\nThank you for using Keshav's Genetic Analysis Program. Goodbye!\n";
		goto program_end;


	}else if(menu_option == '1'){ //Enter Sequence
		enter_sequence();
		
		cout << "Done.";
		getline(cin, dummy);
		goto menu_start;


	}else if(menu_option == '2'){ //Print Sequence/Complement
		if(output_type == '0'){
			get_output_type();
		}else{}
		cout << "\nSequence:\n";
		print_nt(sequence);

		cout << "Starting sequence done.\n\n";
		getline(cin, dummy);
		cout << "  *   *   *   *   *   *   *\n"
			 << "    *   *   *   *   *   *  \n"
			 << "  *   *   *   *   *   *   *\n";

		cout << "\nComplement:\n";
		print_nt(complement);
		cout << endl << "Done.";
		getline(cin, dummy);

		goto menu_start;


	}else if(menu_option == '3'){ //Find ORFs
		int min_ORF_size = 0;
		cout << "What is the minimum size ORF you want (in nucleotides)?\n"
			<<	"Enter 0 (zero) to print *all* ORFs. (WARNING: This may create very large "
			<<	"outputs for even moderately sized fragments.)\n";
		cin >> min_ORF_size;

		//Protecting against invalid inputs for min_ORF_size.
		//If it is alphanumeric but not alphabetic, then the program will proceed as normal.
		//Otherwise, it is assumed that the input is invalid, and print all ORFs.
		if(isalnum(min_ORF_size) != 0 && isalpha(min_ORF_size) == 0){
		}else{
			cout << "Invalid input. By default, all ORFs will be printed.\n";
			min_ORF_size = 0;
		}

		cout << "\n\n\nOriginal Strand:\n";
		ORF_finder(sequence, min_ORF_size);
		
		cout << "*   *   *   *   *   *   *\n"
			 << "  *   *   *   *   *   *  \n"
			 << "*   *   *   *   *   *   *\n";

		cout << "Original Strand Complete.";
		getline(cin, dummy);
		
		cout << endl << endl << endl << "Complementary Strand:\n";
		ORF_finder(complement, min_ORF_size);

		cout << "Done.";
		getline(cin, dummy);

		goto menu_start;


	}else if(menu_option == '4'){ //Get AA seq'n
		if(output_type == '0'){
			get_output_type();
		}else{}
		cout << "\nOriginal strand:\n";
		print_polypeptide(sequence);

		getline(cin, dummy);

		cout << "\n\nComplementary strand:\n";
		print_polypeptide(complement);

		cout << "\nDone.";
		getline(cin, dummy);

		goto menu_start;

	}else if(menu_option == '5'){ //Select Output Type
		get_output_type();
		cout << "\nDone.";
		getline(cin, dummy);
		goto menu_start;


	}else{ //Default
		cout << "Invalid Entry.";
		cin >> dummy;
		goto menu_start;
	}

	
program_end:
	return program_success_code;
}
//END FUCNTION MAIN()


//START ACCESSORY FUNCTIONS


//ORF_FINDER
void ORF_finder(string sequence, int minimum){
	int ORF_start, ORF_end;
	int ORF_length = 0;
	int sequence_length = sequence.size();
	string prot_seq;
	bool is_min;


	for(int frame = 0; frame < 3; frame++){

		for(int index = frame; index < sequence_length; index += 3){
			if(sequence[index] == 'A'){
				if(sequence[index + 1] == 'T'){
					if(sequence[index + 2] == 'G'){
						//int start_site = &sequence[index];
						ORF_start = index + 1;
						prot_seq = ORF_checker(sequence, ORF_start,
							ORF_end, is_min, minimum);
						ORF_length = ORF_end - ORF_start + 1;
						is_min = (ORF_length >= minimum) ? true : false;
						if(is_min){ //CUT: [ORF_end != ORF_start){]
							cout << "\n\n Corresponding DNA sequence is " << ORF_length 
								<< "nt long.\n\n";
							print_nt(sequence.substr(index, ORF_length), "NONE");
							cout << "  *   *   *   *   *   *   *\n"
								 << "    *   *   *   *   *   *  \n"
								 << "  *   *   *   *   *   *   *\n";
						}else {}

						index = ORF_end; //the next ORF must not overlap with this one
					}else {}
				}else {}
			}else{}
		}
	}
}
//END ORF_FINDER


/*
 *Given a sequence and the position of the first nt of an ATG codon, ORF_checker finds
 * the STOP codon in frame with the start site. It then prints the DNA sequence 
 * of the ORF. If there is no STOP codon, the program warns the user and offers to print
 * the sequence anyway. It then returns the index value of the last nt of the STOP codon.
 * If there was no STOP codon, it returns the start site so that the program can find 
 * other ORFs.
 */
//ORF_CHECKER
string ORF_checker(string seqn, int start_site, int& end_site, bool& is_min, int min){
	
	string next_amino = "NONE";
	string potential_protein;
	int length = seqn.size();
	int codon_pos;


	for(codon_pos = (start_site - 1); codon_pos < length; codon_pos += 3){
		next_amino = get_fancy_amino(seqn.substr(codon_pos, 3));
		potential_protein += next_amino;
		if(next_amino == peptide[20]){
			break; //here, codon_pos == position of first nt of STOP codon - 1
		}else {
			potential_protein += " - ";
		}
	}

	end_site = codon_pos + 3; //end_site == position of last nt of STOP codon

	if((end_site - start_site) < min){
		is_min = false;
		goto too_small;
	} else {}
	if(next_amino == peptide[20]){ //if last "amino" was STOP codon
		cout << '(' << start_site << ")\n";
		cout << potential_protein;
		cout << endl << '(' << end_site << ")\n\n";

	}else {	//need to warn user if no STOP codon was found (not a full frame)

		cout << "Grr, Argh! Your sequence contains a start codon that does not have a "
			<< "STOP codon at the end. This is probably not an ORF, unless this "
			<< "fragment is incomplete.\n";
		
		char choice;
		cout << "Would you like to print the frame anyway? Y/N\n";
		cin >> choice;
		if(choice == 'y' || choice == 'Y' || choice == '1'){
			cout << '(' << start_site << ")\n";
			cout << potential_protein;
			cout << endl << '(' << end_site << ")\n";
		}else {}
	}
	
too_small:

	return potential_protein;
}
//END ORF_CHECKER




/*
 * Given a ssDNA sequence, this function will return a complementary sequence.
 *  The return will also be 5' --> 3' (for easier use in the other functions).
 */
//DNA_COMPLEMENT
string DNA_complement(string template_strand){
	string comp;
	int index;
	
	for(index = (length - 1); index >= 0; index--){ //stores the reverse of 
		comp += template_strand[index];				//template_strand to comp
	}

	for(index = 0; index < length; index++){ //performs complementary substitution
		if(comp[index] == 'A'){
			comp[index] = 'T';			// A --> T
		}else if(comp[index] == 'C'){
			comp[index] = 'G';			// C --> G
		}else if(comp[index] == 'G'){
			comp[index] = 'C';			// G --> C
		}else if(comp[index] == 'T'){
			comp[index] = 'A';			// T --> A
		}else {
			comp[index] = '?';
		}
	}
	return comp;

}
//END DNA_COMPLEMENT


/* NOT NEEDED FOR NOW? (with upgrade of simple/fancy get_aminos)
string one_letter_conversion(string AA_seqn){
	int length = AA_seqn.size();
	for(int index = 0; index < length;  ){ //increment is only done if nothing is removed
		if(AA_seqn[index] == ' ' || AA_seqn[index] == '-'){
			AA_seqn.erase(index, 1);
		}else {
			index++;
		}
	}
	length = AA_seqn.size();
	string pep, new_AA_sequence;
	int pep_index;
	for(int index = 0; index < length; index += 3){
		pep = AA_seqn.substr(index, index + 2);
		for(pep_index = 0; pep_index < size_of_arrays; pep_index++){
			if(pep == peptide[pep_index]){		//searching for pep
				break;				//found pep, pep_index has index of AA
			}else{}
		}
		new_AA_sequence += one_letter[pep_index];
	}

	
	return new_AA_sequence;
}
*/


void print_pretty_DNA(string seqn){
	int length = seqn.size();
	int nt_count = 0;
	while(nt_count < length){
		cout << (nt_count + 1) << "\t";
		for(int block = 1; block <= 6; block++){
			if(seqn.substr(nt_count).size() >= 10){
				cout << seqn.substr(nt_count, 10);
				cout << " ";
				nt_count += 10;
			}else{
				cout << seqn.substr(nt_count);
				nt_count += seqn.substr(nt_count).size();
				break;
			}
		}
		cout << endl;
	}
}



/*
 * As you might expect, the print_nt function accepts a DNA sequence and prints
 *  prints it, or calls the pretty print function if needed.
 */
//PRINT_NT
void print_nt(string seqn, string comp){
	if(output_type == '1'){
		cout << "Your cDNA sequence:\n\n5\'\n";
		print_pretty_DNA(seqn);
		cout << "3\'\n";
		cout << endl << endl;
		if(comp != "NONE"){
			cout << "Complementary strand:\n\n5\'\n";
			print_pretty_DNA(comp);
			cout << "3\'\n";
		}
	}else{
		cout << "Your cDNA sequence:\n\n5\'\n" << seqn;
		cout << endl << endl;
		cout << "3\'\n";
		if(comp != "NONE"){
			cout << "Complementary strand:\n\n5\'\n" << comp;
			cout << "\n3\'\n";
		}
	}
}
//END PRINT_NT


/*
 * As you might expect, the print_polypeptide function accepts a string and prints
 *  the corresponding amino acid translation.
 */
//PRINT_POLYPEPTIDE
void print_polypeptide(string seqn){
	cout << "Amino Acid Sequence:\n";
	int length = seqn.size();
	string next_codon;
	string AA_seqn = "";
	int index = 0;

	//Building the AA_seqn string, including divisions into AA fragments.
start_new:
	AA_seqn += (index + 1);
	AA_seqn += "\n";
	while(index < length){
		
		if(output_type == '1'){
			AA_seqn += " - ";
			next_codon = get_fancy_amino(seqn.substr(index, 3)); //local to for loop
			AA_seqn += next_codon;

		}else if(output_type == '2'){
			next_codon = get_simple_amino(seqn.substr(index, 3)); //local to for loop
			AA_seqn += next_codon;
		}

		if(next_codon == peptide[20] || next_codon == one_letter[20]){
			AA_seqn += "\n";
			AA_seqn += (index + 3);
			AA_seqn += "\n\n";
			index += 3;
			goto start_new;
		}else {}
	
		index += 3;
	}

	cout << AA_seqn;
}
//END PRINT_POLYPEPTIDE



void enter_sequence(){
	string new_line;
	string dump;
	getline(cin, dump);
	cout << "Please enter a DNA sequence (without header):\n";

	do{
		getline(cin, new_line);
		sequence += new_line;
	}while(new_line != "\n" && new_line != "\r" && new_line != "");

	length = int(sequence.size()); //need length for loop

	for(int index_format = 0; index_format < length;  ){ //increment ONLY if input is valid
		if(isalpha(sequence[index_format])){
			sequence[index_format] = toupper(sequence[index_format]);
			index_format++;
		}else{
			sequence.erase(index_format, 0);
		}
	}

	length = sequence.size(); //need to find the length of the sequence without spaces

	complement = DNA_complement(sequence);
}

void get_output_type(){
	do{
		cout << "What kind of format do you want for output? (1 = Fancy, 2 = Practical)\n";
		cin >> output_type;
	}while(output_type != '1' && output_type != '2');
}


//GET_SIMPLE_AMINO
string get_simple_amino(string codon){
	string amino;
	if(codon[0] == 'T'){ //Tnn

		if		(codon[1] == 'T'){ //TTn
			if(codon[2] == 'T' || codon[2] == 'C'){ 
				amino = one_letter[13]; //Phenylalanine
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[10]; //Leucine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'C'){ //TCn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[15]; //Serine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'A'){ //TAn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = one_letter[18]; //Tyrosine
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[20]; //STOP
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'G'){ //TGn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = one_letter[4]; //Cysteine
			}else if(codon[2] == 'A'){
				amino = one_letter[20]; //STOP
			}else if(codon[2] == 'G'){
				amino = one_letter[17]; //Tryptophan
			}else {
				amino = one_letter[21]; //ERROR
			}
		}else { //T??
			amino = one_letter[21]; //ERROR
		}


	}else if(codon[0] == 'C'){ //Cnn
		
		if(codon[1] == 'T'){ //CTn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[10]; //Leucine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'C'){ //CCn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[14]; //Proline
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'A'){ //CAn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = one_letter[8]; //Histidine
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[6]; //Glutamine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'G'){ //CGn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[1]; //Arginine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else { //C??
			amino = one_letter[21]; //ERROR
		}


	}else if(codon[0] == 'A'){ //Ann

		if(codon[1] == 'T'){ //ATn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A'){ 
				amino = one_letter[9]; //Isoleucine
			}else if(codon[2] == 'G'){
				amino = one_letter[12]; //Methionine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'C'){ //ACn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[16]; //Threonine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'A'){ //AAn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = one_letter[2]; //Asparagine
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[11]; //Lysine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'G'){ //AGn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = one_letter[15]; //Serine
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[1]; //Arginine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else { //A??
			amino = one_letter[21]; //ERROR
		}

	}else if(codon[0] == 'G'){ //Gnn
		
		if(codon[1] == 'T'){ //GTn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[19]; //Valine
			}else {
				amino = one_letter[21]; //ERROR GT?
			}
			
		}else if(codon[1] == 'C'){ //GCn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[0]; //Alanine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'A'){ //GAn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = one_letter[3]; //Aspartic Acid
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[5]; //Glutamic Acid
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else if(codon[1] == 'G'){ //GGn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = one_letter[7]; //Glycine
			}else {
				amino = one_letter[21]; //ERROR
			}

		}else { //G??
			amino = one_letter[21]; //ERROR
		}

	}else { //???
		amino = one_letter[21]; //ERROR
	}

	return amino;

}
//END GET_SIMPLE_AMINO



//GET_FANCY_AMINO
string get_fancy_amino(string codon){
	string amino;
	if(codon[0] == 'T'){ //Tnn

		if		(codon[1] == 'T'){ //TTn
			if(codon[2] == 'T' || codon[2] == 'C'){ 
				amino = peptide[13]; //Phenylalanine
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[10]; //Leucine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'C'){ //TCn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[15]; //Serine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'A'){ //TAn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = peptide[18]; //Tyrosine
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[20]; //STOP
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'G'){ //TGn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = peptide[4]; //Cysteine
			}else if(codon[2] == 'A'){
				amino = peptide[20]; //STOP
			}else if(codon[2] == 'G'){
				amino = peptide[17]; //Tryptophan
			}else {
				amino = peptide[21]; //ERROR
			}
		}else { //T??
			amino = peptide[21]; //ERROR
		}


	}else if(codon[0] == 'C'){ //Cnn
		
		if(codon[1] == 'T'){ //CTn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[10]; //Leucine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'C'){ //CCn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[14]; //Proline
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'A'){ //CAn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = peptide[8]; //Histidine
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[6]; //Glutamine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'G'){ //CGn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[1]; //Arginine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else { //C??
			amino = peptide[21]; //ERROR
		}


	}else if(codon[0] == 'A'){ //Ann

		if(codon[1] == 'T'){ //ATn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A'){ 
				amino = peptide[9]; //Isoleucine
			}else if(codon[2] == 'G'){
				amino = peptide[12]; //Methionine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'C'){ //ACn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[16]; //Threonine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'A'){ //AAn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = peptide[2]; //Asparagine
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[11]; //Lysine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'G'){ //AGn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = peptide[15]; //Serine
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[1]; //Arginine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else { //A??
			amino = peptide[21]; //ERROR
		}

	}else if(codon[0] == 'G'){ //Gnn
		
		if(codon[1] == 'T'){ //GTn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[19]; //Valine
			}else {
				amino = peptide[21]; //ERROR GT?
			}
			
		}else if(codon[1] == 'C'){ //GCn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[0]; //Alanine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'A'){ //GAn
			if(codon[2] == 'T' || codon[2] == 'C'){
				amino = peptide[3]; //Aspartic Acid
			}else if(codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[5]; //Glutamic Acid
			}else {
				amino = peptide[21]; //ERROR
			}

		}else if(codon[1] == 'G'){ //GGn
			if(codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A' || codon[2] == 'G'){
				amino = peptide[7]; //Glycine
			}else {
				amino = peptide[21]; //ERROR
			}

		}else { //G??
			amino = peptide[21]; //ERROR
		}

	}else { //???
		amino = peptide[21]; //ERROR
	}

	return amino;

}
//END GET_FANCY_AMINO
