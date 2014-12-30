# Usage
#sh format.sh data/database.txt "C[C@H](O)C(O)=O"
#

java -jar ChemBlast/dist/ChemBlast.jar -searchDB $1 -querySMI $2
