from flask import Flask, request, url_for, jsonify, json, render_template, send_file#make_response
from flask_restful import reqparse, abort, Api, Resource
#from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts, Draw
from StringIO import StringIO
import PAINSrelief
#from PAINSrelief import pains2

app = Flask(__name__)


@app.route('/_PAINS',methods=['POST'])
def PAINS():
    chembenchargs=request.get_json() #Turn the request into a JSON for easy management.
    smiles=chembenchargs['smiles'] #Pull the smile out of the post request.
    print('Request being made for ' + smiles)
    martin = PAINSrelief.pains(smiles) #Returns a JSON object.
    response = martin
    return response #Send the JSON object to the browser.

@app.route('/', methods=['GET'])
def index():
    return app.send_static_file('index.html') #Whne querying the first time, send them our home page.

@app.route('/MCRA/mcra.html')
def MCRA():
    return app.send_static_file('mcra.html')

@app.route('/image/<smarts_or_smile>/<smile>.png')
def image(smile, smarts_or_smile):
    print(smile)
    if(smarts_or_smile == 'smarts'): m = MolFromSmarts(str(smile))
    else: m = MolFromSmiles(str(smile))
    if(m==None):
        #If our inputted smile is poorly formatted we just return a simple molecule	
         m = MolFromSmiles('c1ccccc1')	 
    img = StringIO() #A file buffer to allow us to return the image.
    m.UpdatePropertyCache(strict=False)
    pil = Draw.MolToImage(m,kekulize=False) #Pillow image generated from our mol.
    pil.save(img, 'PNG', quality=100) #Put the pillow image in our StringIO
    img.seek(0) #Reset our buffer to 0.
    return send_file(img, mimetype='image/png')

@app.route('/highlighted-image/<smarts_or_smile>/<inds>_<smile>.png')
def highlighted_image(smile, smarts_or_smile, inds):
    ind_array = [int(x) for x in inds.split("-")]
    smile = str(smile).replace('FOWARDSLASH', '/')
    smile = str(smile).replace('UNDERSCORE', '_')
    if(smarts_or_smile == 'smarts'): m = MolFromSmarts(str(smile))
    else: m = MolFromSmiles(str(smile))
    if(m==None):
        #If our inputted smile is poorly formatted we just return a simple molecule	
	     m = MolFromSmiles('c1ccccc1') 
    img = StringIO() #A file buffer to allow us to return the image.
    pil = Draw.MolToImage(m, highlightAtoms=ind_array) #Pillow image generated from our mol.
    pil.save(img, 'PNG', quality=100) #Put the pillow image in our StringIO
    img.seek(0) #Reset our buffer to 0.
    return send_file(img, mimetype='image/png')

@app.route('/<path:path>')
def jsme(path):
	print path
	return app.send_static_file(path)

if __name__ == '__main__':
	app.run(host='0.0.0.0',port=80)
