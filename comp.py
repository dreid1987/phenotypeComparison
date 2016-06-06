import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import hypergeom
import datetime
today=datetime.datetime.today()

allOrthologues=False
countOMIM=False
comparePhenotype=False
compareMutationImpact=True
compareKnockoutMorpho=False
compareMorphoCRISPR=False

class mdict(dict):
	def __setitem__(self, key, value):
		 self.setdefault(key, []).append(value)
		 
		 
def loadSpeciesOrthos(species):
	
	homToSpec=mdict()
	specToHom=dict()
	
	for line in open('database/HOM_AllOrganism.rpt'):
		line=line.split('\t')
		if len(line)>1:
			if line[1]==species:
				hom=line[0]
				gene=line[3]
				homToSpec[hom]=gene
				specToHom[gene]=hom
	return homToSpec,specToHom

if allOrthologues:
	mou=[]
	human=[]
	fish=[]
	for line in open('database/HOM_AllOrganism.rpt'):
		line=line.split('\t')
		if len(line)>1:
			if line[1]=='mouse, laboratory':
				mou.append(line[0])
			if line[1]=='human':
				human.append(line[0])
			if line[1]=='zebrafish':
				fish.append(line[0])
	print '      \tHum\tMou\tFish'
	print 'All   ' + '\t' + str(len(human)) + '\t' + str(len(mou)) + '\t' + str(len(fish))
	
	for species in ['mou','human','fish']:
		mouC=0
		humanC=0
		fishC=0
		for gene in vars()[species]:
			if gene in human:
				humanC+=1
			if gene in mou:
				mouC+=1
			if gene in fish:
				fishC+=1
		
		print 'In ' + species +'\t' + str(humanC) + '\t' + str(mouC) + '\t' + str(fishC)
homToMouse,mouseToHom=loadSpeciesOrthos('mouse, laboratory')
homToHuman,humanToHom=loadSpeciesOrthos('human')
homToFish,fishToHom=loadSpeciesOrthos('zebrafish')
if countOMIM:
	
	"""
	lens=[]
	for gene in homToHuman.keys():
		lens.append(len(homToHuman[gene]))
	
	for i in range(5):
		
		print i, lens.count(i)
	"""
	
	
	fishModels=[]
	for line in open('database/phenoGeneCleanData_fish.txt'):
		line=line.split('\t')
		
		gene=line[1]
		
		try:
			homID=fishToHom[gene]
			fishModels.append(homID)
		except KeyError: pass
	mouseModels=[]
	for line in open('database/MGI_PhenoGenoMP.rpt'):
		line=line.split('<')
		try:
			gene=line[0]
			homID=mouseToHom[gene]
			mouseModels.append(homID)
		except KeyError: pass
	
	humanDiseaseGenes=[]
	for line in open('database/genemap.txt'):
		if line[0]!='#':
		
			line=line.split('\t')
			
			genes = line[5].split(',')
			for gene in genes:
				try:
					homID=humanToHom[gene]
					humanDiseaseGenes.append(homID)
				except KeyError: pass
		
	humanDiseaseGenes=set(humanDiseaseGenes)
	
	inFish=0
	inMouse=0
	inHuman=0
	inAll=0
	
	for humanGene in humanDiseaseGenes:
		inHuman+=1
		fish=False
		mouse=False
		if humanGene in fishModels:
			fish=True
			inFish+=1
		if humanGene in mouseModels:
			mouse=True
			inMouse+=1
		if fish and mouse:
			inAll+=1
	print inHuman,inMouse,inFish,inAll

if comparePhenotype:	
	def getPubYearDiffsFish(fishYearDiffsx,fishYearDiffsy,gene,humanYear,zfishAnatomy,zfishPhenotype):
		pubs=geneToPubs[gene]
		pubs=set(pubs)
		for pub in pubs:
			try:
				yearDiff=pubYear[pub]-humanYear
				foundPhenotype=0
				for phenotype in pubToPhenotypes[pub]:
					
					p=0
					while p<len(zfishAnatomy):
						phenotypeOK=False
						if phenotype[1]==zfishPhenotype[p] or zfishPhenotype[p]=='':
							phenotypeOK=True
						if phenotype[0]==zfishAnatomy[p] and phenotypeOK:
							foundPhenotype=1
						p+=1
				fishYearDiffsx.append(yearDiff)
				fishYearDiffsy.append(foundPhenotype)
			except KeyError: pass
		return 	fishYearDiffsx,fishYearDiffsy
	#Get fish publication data
	pubYear=dict()
	for line in open('database/zfinpubs.txt'):
		if len(line)>4 and line[0] != '#':
			
			line=line.split('\t')
			year=int(line[5])
			pubID=line[0]
			pubYear[pubID]=year
	
	
	#Get fish microcephaly frequency:
	hasMicrocephaly=0;all=0.
	lastGene=''
	fishMicrocephaly=[]
	fishAllGenes=[]
	fishPhenotypes=mdict()
	
	geneToPubs=mdict()		
	pubToPhenotypes=mdict()
	
	zfishAnatomy=['ZFA:0000100']
	zfishPhenotype=[''] #zero-length string to allow for any anatomy annotation
	mousePhenotypesReference=['MP:0000851','MP:0000849','MP:0000850','MP:0000852','MP:0008886']
	humanPhenotype=['cerebell']
	humanPhenotypeName=''
	try:
		dirTouse='plots/comparePhenotype/'
		for phen in humanPhenotype:
			dirTouse=dirTouse+phen
			humanPhenotypeName =  humanPhenotypeName + phen
		os.mkdir(dirTouse)
	except OSError: pass
	
	for line in open('database/phenoGeneCleanData_fish.txt'):
		
		if line[0]!='#':
			line=line.split('\t')
			gene=line[1]
			pub=line[23]
			pheno=[line[7],line[9]]
			try:
				hom=fishToHom[gene]
				pubToPhenotypes[pub]=pheno
				p=0
				while p<len(zfishPhenotype):
					phenotypeOK=False
					if line[9]==zfishPhenotype[p] or zfishPhenotype[p]=='':
						phenotypeOK=True
					if phenotypeOK and line[7]==zfishAnatomy[p]:
						fishMicrocephaly.append(hom)
						break
					p+=1
				
				phenotype= line[8] + ' ' + line[10] + ' ' + line[11] + ' ' + pub
				fishPhenotypes[hom]=phenotype
				fishAllGenes.append(hom)
				geneToPubs[hom]=pub
				
				
			except KeyError: pass		
	numFishMicrocephaly=len(set(fishMicrocephaly))			
	numFishGenes=len(set(fishAllGenes))		
	print 'Fish: ', str(100.*numFishMicrocephaly/numFishGenes), numFishMicrocephaly, numFishGenes
	
	
	#Get mouse microcephaly frequency:
	mousePhenoDefinition=dict()
	for line in open('database/VOC_MammalianPhenotype.rpt'):
		if len(line)>2 and line[0] != '#':
			line=line.split('\t')
			mousePhenoDefinition[line[0]]=line[1]
	
	hasMicrocephaly=0;all=0.
	usedGenes=[]
	mouseMicrocephaly=[]
	mouseAllGenes=[]
	mousePhenotypes=mdict()
	for line in open('database/MGI_PhenoGenoMP.rpt'):
		if line[0]!='#':
			line=line.strip('\n').split('\t')
			gene=line[1].split('<')[0]
			try:
				hom=mouseToHom[gene]
				mousePhenotypes[hom] = mousePhenoDefinition[line[3]] + ' ' + line[3] + ' PMID:' + line[4]
				for mousePhenotype in mousePhenotypesReference:
					if line[3]== mousePhenotype:
						mouseMicrocephaly.append(hom)
						break
				
				
				mouseAllGenes.append(hom)
			except KeyError: pass	
				
	numMouseMicrocephaly=len(set(mouseMicrocephaly))			
	numMouseGenes=len(set(mouseAllGenes))	
	
	print 'Mouse: ', str(100.*numMouseMicrocephaly/numMouseGenes), numMouseMicrocephaly, numMouseGenes
	
	geneInMouse=0
	mouseHasMicrocephaly=0
	
	geneInFish=0
	fishHasMicrocephaly=0
	humanMicrocephalyGenes=[]
	usedGenes=[]
	fishYearDiffsx=[]
	fishYearDiffsy=[]
	
	
	
	
	print 'Gene______'+ '\t' + 'Fish' + '\t' + 'Mouse' + '\t' + 'humanPhenotype'
	print '----------------------------------------------------------------------'
	
	toWrite='\n#Gene      '+ '\t' + 'Fish' + '\t' + 'Mouse' + '\t' + 'humanPhenotype' + '\n#----------------------------------------------------------------------'
	for line in open('database/genemap.txt'):
		if line[0]!='#':
		
			line=line.split('\t')
			
			phenotype=line[11].lower()
			try:
				phenoFound=False
				for refPheno in humanPhenotype:
					
				
					if phenotype.find(refPheno)>-1:
						phenoFound=True
				if phenoFound:
					
					genes=line[5].split(',')
					humanYear=int(line[3])
					if humanYear>16:
						humanYear+=1900
					else:
						humanYear+=2000
					for gene in genes:
						humanGene=humanToHom[gene]
						if humanGene not in usedGenes:
							
							usedGenes.append(humanGene)
							humanMicrocephalyGenes.append(humanToHom[gene])
							fishFound="x";mouseFound="x"
							
							if humanGene in fishAllGenes:
								geneInFish+=1
								if humanGene in  fishMicrocephaly:
									fishHasMicrocephaly+=1
									fishFound=True
									fishYearDiffsx,fishYearDiffsy=getPubYearDiffsFish(fishYearDiffsx,fishYearDiffsy,humanGene,humanYear,zfishAnatomy,zfishPhenotype)
								else:
									fishFound=False
							if humanGene in mouseAllGenes:
								geneInMouse+=1
								if humanGene in  mouseMicrocephaly:
									mouseHasMicrocephaly+=1
									mouseFound=True
									#print 'mouse ', gene
							
								else:
									mouseFound=False
							
							print gene + ' '*(10-len(gene)) + '\t' + str(fishFound) + '\t' + str(mouseFound)+'\t' + phenotype 
							toWrite= toWrite+'\n' + gene+ ' '*(10-len(gene)) +  '\t' + str(fishFound) + '\t' + str(mouseFound)+'\t' + phenotype.replace('?','')
							if fishFound in [True,False]:
								ptypes=set(fishPhenotypes[humanGene])
								for phenotype in ptypes:
									print '\t' + phenotype
									toWrite = toWrite + '\n\t\tFish:' + phenotype
							if mouseFound in [True,False]:
								ptypes=set(mousePhenotypes[humanGene])
								for phenotype in ptypes:
									print '\t' + phenotype
									toWrite = toWrite + '\n\t\t\tMouse:' + phenotype
							toWrite=toWrite + '\n'
			except KeyError: pass
	
	write=open(dirTouse + '/' + humanPhenotypeName+'Table.txt','w')
	
	write.writelines('#Human phenotype: ' + humanPhenotypeName)
	write.writelines('\n#\tAnalysis run ' + str(today.year) + '-' + str(today.month) + '-' + str(today.day))
	
	
	write.writelines('\n#\t' +  str(len(set(humanMicrocephalyGenes))) + ' genes of ' + str(len(humanToHom)) + ' give rise to this phenotype')
	
	
	write.writelines('\n#Zebrafish phenotype: ' + zfishAnatomy[0]+ ' ' + zfishPhenotype[0])
	p=1
	while p<len(zfishAnatomy):
		 write.writelines(',' + zfishAnatomy[p]+ ' ' + zfishPhenotype[p])
		 p+=1
	write.writelines('\n#\t' + str(numFishMicrocephaly) + ' of ' + str(numFishGenes) + ' fish genes match phenotype (' + str(100.*numFishMicrocephaly/numFishGenes) + '%)')
	write.writelines('\n#\t' + str(fishHasMicrocephaly) + ' of ' + str(geneInFish) + ' fish with human genes match phenotype (' + str(100.*fishHasMicrocephaly/geneInFish) + '%)') 
	fishPval=hypergeom.pmf(fishHasMicrocephaly,numFishGenes,numFishMicrocephaly,geneInFish)
	write.writelines('\n#\tp-value = ' + str(fishPval) + ' by hypergeometic test') 
	
	write.writelines('\n#Mouse phenotype: ' + mousePhenotypesReference[0])
	for pheno in mousePhenotypesReference[1:]:
		write.writelines(',' + pheno)
	write.writelines('\n#\t' + str(numMouseMicrocephaly) + ' of ' + str(numMouseGenes) + ' mice genes match phenotype (' + str(100.*numMouseMicrocephaly/numMouseGenes) + '%)')
	write.writelines('\n#\t' + str(mouseHasMicrocephaly) + ' of ' + str(geneInMouse) + ' mice with human genes match phenotype (' + str(100.*mouseHasMicrocephaly/geneInMouse) + '%)') 
	mousePval=hypergeom.pmf(mouseHasMicrocephaly,numMouseGenes,numMouseMicrocephaly,geneInMouse)
	write.writelines('\n#\tp-value = ' + str(mousePval) + ' by hypergeometic test') 
	
	write.writelines(toWrite)
	
	
	write.close()
	if fishHasMicrocephaly>0:
		#plt.plot(fishYearDiffsx,fishYearDiffsy,'o')
		import pandas as pd
		import statsmodels.api as sm
		import patsy
		X=np.array(fishYearDiffsx)
		y=np.array(fishYearDiffsy)
		print len(X),len(y)
	
		from sklearn import linear_model
		X = X[:, np.newaxis]

		# run the classifier
		clf = linear_model.LogisticRegression(C=1e5)
		clf.fit(X, y)

		# and plot the result
		plt.figure(1, figsize=(4, 3))
		plt.clf()
		plt.scatter(X.ravel(), y, color='black', zorder=20)
		X_test = np.linspace(min(X), max(X), 300)


		def model(x):
		    return 1 / (1 + np.exp(-x))
		loss = model(X_test * clf.coef_ + clf.intercept_).ravel()
		plt.plot(X_test, loss, color='blue', linewidth=3)

		ols = linear_model.LinearRegression()
		ols.fit(X, y)
	
		plt.plot(X_test, ols.coef_ * X_test + ols.intercept_, linewidth=1)
		plt.axhline(.5, color='.5')

	
		plt.ylim(-.25, 1.25)
		plt.xlim(min(X), max(X))
		plt.xlabel('Year relative to human discovery')
		plt.yticks(())
		plt.ylabel('Phenotype identified in fish?')
		plt.savefig(dirTouse + '/'+humanPhenotypeName+'pubYear.svg')
		plt.savefig(dirTouse + '/'+humanPhenotypeName+'pubYear.png')
		plt.clf()
def overlappingSet(list1,list2):
	list1=set(list1)
	list2=set(list2)
	
	list1Only=list1.difference(list2)
	list2Only=list2.difference(list1)
	both=list1.intersection(list2)
	
	return list1Only, list2Only, both
def getPvalHypergeom(go,tot1s,totGenes):
	
	hit=0;tot=0
	for gene in GOs[go]:
	
		try:
			hit+=scores[gene]
			tot+=1	
		except KeyError:
			fail=True
	
	stat=hypergeom.sf(int(hit),totGenes,tot1s,tot,loc=0)
	
if compareMutationImpact:	
	deletionsForGene=mdict()
	variantsForGene=mdict()
	mutationType=dict()
	geneName=dict()
	for line in open('database/genotype_features.txt'):
		if len(line)>5:
			line=line.split('\t')
			gene = line[9]
			geneName[gene]=line[8]
			if len(gene)>0:
				genotype=line[0]
				#if line[6]=='deletion' or line[8]=='deficiency':
				if line[8]=='deficiency':
					deletionsForGene[gene]=genotype
					mutationType[genotype]='Deletion'
				if line[7]=='Point Mutation' or line[6]=='sequence variant':
					variantsForGene[gene]=genotype
					mutationType[genotype]='Point mutation'
				
				
	for line in open('database/genotype_features_missing_markers.txt'):
		if len(line)>5:
			line=line.split('\t')
			gene=line[4]
			genotype=line[0]
			deletionsForGene[gene]=genotype
			mutationType[genotype]='Deletion'			

	
	fishForGenotype=mdict()#Exclude morpholino-treated fish
	fishGenotype=dict()
	for line in open('database/fish_components_fish.txt'):
		line=line.split('\t')
		if line[4][:7]=='ZDB-ALT' and line[1].find('+MO')==-1 and line[1].find('+ MO')==-1:
			fishForGenotype[line[10]]=line[0]
			fishGenotype[line[0]]=line[1]
	
	phenotypeForFish=mdict()
	phenotypeList=[]
	for line in open('database/phenotype_fish.txt'):
		line=line.split('\t')
		
		phenotype=line[11]
		for spot in line[13:18]:
			if len(spot)>0:
				phenotype = phenotype + ' ' + spot
		phenotypeForFish[line[0]]=phenotype
		phenotypeList.append(phenotype)
	numPhenotypes=len(set(phenotypeList))
	
	genesToUse=[]
	for gene in deletionsForGene.keys():
		foundDeletionFish=False
		for deletion in deletionsForGene[gene]:
			if deletion in fishForGenotype.keys():
				foundDeletionFish=True
		
		if gene in variantsForGene.keys():
			foundVariantFish=False
			for variant in variantsForGene[gene]:
				if variant in fishForGenotype.keys():
					foundVariantFish=True
			
			if foundDeletionFish and foundVariantFish:
				genesToUse.append(gene)

		
	
	print len(genesToUse)
	write=open('plots/pointVsKO.txt','w')
	bothPercent=[];delOnlyPercent=[];varOnlyPercent=[]
	highPs=0;lowPs=0
	for gene in genesToUse:
		
	
		
		variantPhenotypes=[]
		deletionPhenotypes=[]
		writeFish=''	
		for genotype in set(variantsForGene[gene] + deletionsForGene[gene]):
			try:
				for fish in set(fishForGenotype[genotype]):
				
					writeFish=writeFish + '\n\t' +mutationType[genotype] + ' - ' + fish + ': ' + fishGenotype[fish] + '\n\t\t'
				
					try:
						for phenotype in phenotypeForFish[fish]:
							writeFish=writeFish  + phenotype + ','
							if mutationType[genotype] == 'Point mutation':
								variantPhenotypes.append(phenotype)
							if mutationType[genotype] == 'Deletion':
								deletionPhenotypes.append(phenotype)
					except KeyError:
						writeFish=writeFish + 'No annotated phenotype'
			except KeyError: pass
		variantOnly,deletionOnly,both=overlappingSet(variantPhenotypes,deletionPhenotypes)		
		total=float(np.sum([len(variantOnly),len(deletionOnly),len(both)]))
		
		if len(variantPhenotypes)>0 and len(deletionPhenotypes)>0:
			overlapPvalue=hypergeom.pmf(len(both),numPhenotypes,len(both)+len(variantOnly),len(both)+len(deletionOnly))
			if overlapPvalue<0.05:
				lowPs+=1
			else:
				highPs+=1
			write.writelines('\n\n' + geneName[gene] + '\t' + gene + '\tP=' + str(overlapPvalue))
			write.writelines('\n\t' + str(len(deletionOnly)) + ' phenotypes in deletion only:')
			for pheno in deletionOnly:
				write.writelines(', ' + pheno)
			write.writelines('\n\t' + str(len(variantOnly)) + ' phenotypes in variants only:')
			for pheno in variantOnly:
				write.writelines(', ' + pheno)
			write.writelines('\n\t' + str(len(both)) + ' phenotypes in both:')
			for pheno in both:
				write.writelines(', ' + pheno)	
			write.writelines('\n' + writeFish)	
			
			
			
			bothPercent.append(len(both)/total)
			delOnlyPercent.append(len(deletionOnly)/total)
			varOnlyPercent.append(len(variantOnly)/total)
		
	write.close()			
			
	print np.mean(delOnlyPercent)*100,np.mean(bothPercent)*100,np.mean(varOnlyPercent)*100
	print lowPs,highPs
if compareKnockoutMorpho:
	deletionsForGene=mdict()
	
	mutationType=dict()
	geneName=dict()
	for line in open('database/genotype_features.txt'):
		if len(line)>5:
			line=line.split('\t')
			gene = line[9]
			geneName[gene]=line[8]
			if len(gene)>0:
				genotype=line[0]
				#if line[6]=='deletion' or line[8]=='deficiency':
				if line[8]=='deficiency':
					deletionsForGene[gene]=genotype
					mutationType[genotype]='Deletion'
				
				
				
	for line in open('database/genotype_features_missing_markers.txt'):
		if len(line)>5:
			line=line.split('\t')
			gene=line[4]
			genotype=line[0]
			deletionsForGene[gene]=genotype
			mutationType[genotype]='Deletion'			

	
	fishForGenotype=mdict()#Exclude morpholino-treated fish
	fishForGeneMorpho=mdict()
	fishGenotype=dict()
	morphoToGene=dict()
	for line in open('database/Morpholinos.txt'):	
		line=line.split('\t')
		morphoToGene[line[3]]=line[0]
		
		
	for line in open('database/fish_components_fish.txt'):
		line=line.split('\t')
		fishGenotype[line[0]]=line[1]
		if line[4][:7]=='ZDB-ALT' and line[1].find('+MO')==-1 and line[1].find('+ MO')==-1:
			fishForGenotype[line[10]]=line[0]
			
		if line[1].find('+MO')>-1 or line[1].find('+ MO')>-1:
			try:
				gene=morphoToGene[line[4]]
				fishForGeneMorpho[gene]=line[0]
			except KeyError: pass
	
	phenotypeForFish=mdict()
	phenotypeList=[]
	for line in open('database/phenotype_fish.txt'):
		line=line.split('\t')
		
		phenotype=line[11]
		for spot in line[13:18]:
			if len(spot)>0:
				phenotype = phenotype + ' ' + spot
		phenotypeForFish[line[0]]=phenotype
		phenotypeList.append(phenotype)
	numPhenotypes=len(set(phenotypeList))
	
	genesToUse=[]
	for gene in deletionsForGene.keys():
		for deletion in deletionsForGene[gene]:
			if deletion in fishForGenotype.keys():
				if gene in fishForGeneMorpho.keys():
					genesToUse.append(gene)
		

		
	
	write=open('plots/morphoVsKO.txt','w')
	bothPercent=[];delOnlyPercent=[];morphoOnlyPercent=[]
	highPs=0;lowPs=0
	for gene in set(genesToUse):
		
		deletionFish=[]
		for genotype in deletionsForGene[gene]:
			for fish in set(fishForGenotype[genotype]):
				deletionFish.append(fish)
		
		morphantPhenotypes=[]
		deletionPhenotypes=[]
		writeFish=''	
		for fish in set(deletionFish):

				
			writeFish=writeFish + '\n\tDeletion - ' + fish + ': ' + fishGenotype[fish] + '\n\t\t'
		
			try:
				for phenotype in phenotypeForFish[fish]:
					writeFish=writeFish  + phenotype + ','
					deletionPhenotypes.append(phenotype)
			except KeyError:
				writeFish=writeFish + 'No annotated phenotype'
			except KeyError: pass
			
		for fish in set(fishForGeneMorpho[gene]):
			
			try:
				writeFish=writeFish + '\n\tMorpholino - ' + fish + ': ' + fishGenotype[fish] + '\n\t\t'
		
				try:
					for phenotype in phenotypeForFish[fish]:
						writeFish=writeFish  + phenotype + ','
						morphantPhenotypes.append(phenotype)
				except KeyError:
					writeFish=writeFish + 'No annotated phenotype'
			except KeyError: pass
			
			
		morphoOnly,deletionOnly,both=overlappingSet(morphantPhenotypes,deletionPhenotypes)		
		total=float(np.sum([len(morphoOnly),len(deletionOnly),len(both)]))
		
		if len(deletionPhenotypes)>0 and len(morphantPhenotypes)>0:
			overlapPvalue=hypergeom.pmf(len(both),numPhenotypes,len(both)+len(morphoOnly),len(both)+len(deletionOnly))
			if overlapPvalue<0.05:
				lowPs+=1
			else:
				highPs+=1
			try:
				write.writelines('\n\n' + geneName[gene] + '\t' + gene + '\tP=' + str(overlapPvalue))
			except KeyError: 
				write.writelines('\n\n' + gene + '\tP=' + str(overlapPvalue))
			write.writelines('\n\t' + str(len(deletionOnly)) + ' phenotypes in deletion only:')
			for pheno in deletionOnly:
				write.writelines(', ' + pheno)
			write.writelines('\n\t' + str(len(morphoOnly)) + ' phenotypes in morpholino only:')
			for pheno in morphoOnly:
				write.writelines(', ' + pheno)
			write.writelines('\n\t' + str(len(both)) + ' phenotypes in both:')
			for pheno in both:
				write.writelines(', ' + pheno)	
			write.writelines('\n' + writeFish)	
			
			
			
			bothPercent.append(len(both)/total)
			delOnlyPercent.append(len(deletionOnly)/total)
			morphoOnlyPercent.append(len(morphoOnly)/total)
		
	write.close()			
			
	print np.mean(delOnlyPercent)*100,np.mean(bothPercent)*100,np.mean(morphoOnlyPercent)*100
	print lowPs, highPs

if compareMorphoCRISPR:
		
	geneCRISPRs=mdict()
	geneMorphos=mdict()
	for line in open('database/CRISPR.txt'):
		line=line.split('\t')
		geneCRISPRs[line[0]]=line[3]
	for line in open('database/Morpholinos.txt'):	
		line=line.split('\t')
		geneMorphos[line[0]]=line[3]
	
	
	
	fishForMorpho=mdict()	
	fishForCRISPR=mdict()
	for line in open('database/fish_components_fish.txt'):	
		line=line.split('\t')
		
		if line[4].find('MRPHLNO')>-1:
			fishForMorpho[line[4]]=line[0]
		if line[4].find('CRISPR')>-1:
			fishForCRISPR[line[4]]=line[0]
			
	
	print len(fishForCRISPR.keys())
	
	phenotypeForFish=mdict()
	for line in open('database/phenotype_fish.txt'):
		line=line.split('\t')
		
		phenotype=line[10] + '|' + line[12]
		phenotypeForFish[line[0]]=phenotype
		
	crPercent=[]
	bothPercent=[]
	moPercent=[]
	
	for gene in geneCRISPRs.keys():
		
		if gene in geneMorphos.keys():
			
			try:
				crisprPhenotypes=[]
				crCount=0;moCount=0
				for crispr in geneCRISPRs[gene]:
					
					for fish in fishForCRISPR[crispr]:
						
						for phenotype in phenotypeForFish[fish]:
							
							crisprPhenotypes.append(phenotype)
				moPhenotypes=[]
				for mo in geneMorphos[gene]:
					for fish in fishForMorpho[mo]:
						
						for phenotype in phenotypeForFish[fish]:
							moPhenotypes.append(phenotype)
					
				
				moOnly=0.
				crisprOnly=0.
				both=0.
				for phenotype in moPhenotypes:
					if phenotype in crisprPhenotypes:
						both+=1
					else:
						moOnly+=1
				for phenotype in crisprPhenotypes:
					if phenotype not in moPhenotypes:
						moOnly+=1
			
				s=np.sum([moOnly,crisprOnly,both])
				crPercent.append(crisprOnly/s)
				moPercent.append(moOnly/s)
				bothPercent.append(both/s)
				print moOnly,crisprOnly,both
			except KeyError: pass
	
	print np.median(crPercent)*100,np.median(bothPercent)*100,np.median(moPercent)*100
	
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
