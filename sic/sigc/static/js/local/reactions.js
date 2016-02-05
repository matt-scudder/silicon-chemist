/*
   Takes care of all the stuff related to the test page.
*/
var reactants = [];
var products = [];
var solvent = []; //currently unused
var jsmeApplet = null;

function add_molecule(smiles,mol_type) {
	/*
	   Adds the reactant to the list of reactants, products or solvent (depending on mol_type), 
	   and displays the reactant below the correct column.
	   Right now this will just display SMILES strings. In the future we can use OpenBabel to display images (http://ctr.wikia.com/wiki/Depict_a_compound_as_an_image ) by querying the server with an AJAX request.
	   */
	switch(mol_type){
		case "reactant":
			reactants.push(smiles);
			$("#add_reactant").before("<p>"+smiles+"</p>");
			break;
		case "product":
			products.push(smiles)
			$("#add_product").before("<p>"+smiles+"</p>");
			break;
		case "solvent":
			solvent.push(smiles);
			break;
		default:
			console.error("Bad mol_type passed in to add_molecule: " + mol_type + " for SMILES: " + smiles);
			break;
	}
}

function jsmeOnLoad() {
	jsmeApplet = new JSApplet.JSME("jsme_container","380px","340px"); //figure out better values for the width/height
	$("#add_reactant").click(function() { add_molecule(jsmeApplet.smiles(),"reactant");});
	$("#add_product").click(function() { add_molecule(jsmeApplet.smiles(),"product");});
	//$("#add_solvent").click(function() { add_molecule(jsmeApplet.smiles(),"solvent");});
}

function create_reaction_object(){
	/*
	 * This function generates the JSON for the reaction based on the user's inputted reactants and products.
	 * Separated from send_reaction_to_sic() for modularity.
	 */
       output_object = {
	     reactants:reactants,
	     products:products,
	     solvent: (solvent.length > 0) ? solvent : false
       }; 
       return JSON.stringify(create_reaction_object())
}

function send_reaction_to_sic(){
	$.ajax("/submit_reaction",
			{
				data: create_reaction_object(),
				method:"POST",
				dataType:"text",
				contentType:"application/json",
				success: function(responseData) {
					$("#output_header").show();
					$("#output_container").html(responseData);
				}
			});
}

$("#send_button").click(send_reaction_to_sic);
