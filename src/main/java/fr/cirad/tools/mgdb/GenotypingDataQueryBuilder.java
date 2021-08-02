/*******************************************************************************
 * GIGWA - Genotype Investigator for Genome Wide Analyses
 * Copyright (C) 2016 - 2019, <CIRAD> <IRD>
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License, version 3 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * See <http://www.gnu.org/licenses/agpl.html> for details about GNU General
 * Public License V3.
 *******************************************************************************/
package fr.cirad.tools.mgdb;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import org.apache.commons.math.util.MathUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.bson.types.ObjectId;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

//import com.fasterxml.jackson.databind.ObjectMapper;
import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import com.mongodb.client.MongoCollection;

import fr.cirad.controller.GigwaMethods;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData.VariantRunDataId;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.model.GigwaSearchVariantsRequest;
import fr.cirad.tools.Helper;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class GenotypingDataQueryBuilder.
 */
public class GenotypingDataQueryBuilder implements Iterator<List<BasicDBObject>>
{
	
	/** The Constant LOG. */
	protected static final Logger LOG = Logger.getLogger(GenotypingDataQueryBuilder.class);
	
	/** The mongo template. */
	private MongoTemplate mongoTemplate;


	/** The genotyping project. */
	private GenotypingProject genotypingProject;
	
	/** Whether or not project has effect annotations. */
	private boolean projectHasEffectAnnotations;

	/** The gene names. */
	private String geneNames;
	
	/** The variant effects. */
	private String variantEffects;
	
	/** The selected individuals. */
	private Collection<String>[] selectedIndividuals = new Collection[2];
	
	/** The operator. */
	private String[] operator = new String[2];
	
	/** The percentage of individuals for the "all same" filter. */
	private Integer[] mostSameRatio = new Integer[2];
	
	/** The annotation field thresholds. */
	private HashMap<String, Float>[] annotationFieldThresholds = new HashMap[2];
		
	/** The missing data threshold. */
	private Float[] missingData = new Float[2];
	
	/** The minmaf. */
	private Float[] minmaf = new Float[2];
	
	/** The maxmaf. */
	private Float[] maxmaf = new Float[2];
	
    boolean fDiscriminate = false;
	
	/** The individual index to sample list map. */
	private TreeMap<String /*individual*/, ArrayList<GenotypingSample>>[] individualToSampleListMap = new TreeMap[2];
	
	/** The n total variant count. */
	private long nTotalVariantCount = 0;
	
	private int nNextCallCount = 0;
	
	private int maxAlleleCount = 0;
	
	private List<Integer> filteredGroups;
	
	/** Used for shuffling queried chunks */
	private List<Integer> intervalIndexList;
	
	private List<Map> taggedVariantList;

	private BasicDBList variantQueryDBList;

	private Document groupFields;
	
	private Document projectionFields;

	/** The Constant AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX. */
	static final public String AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX = "_ALL_"; // used to differentiate aggregation query with $and operator 
	
	/** The Constant AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX. */
	static final public String AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX = "_ATLO_";  // used to differentiate find query with $or operator
	
	/** The Constant AGGREGATION_QUERY_NEGATION_SUFFIX. */
	static final public String AGGREGATION_QUERY_NEGATION_SUFFIX = "_NEG_";	// used to indicate that the match operator should be negated in the aggregation query
	
	/** The Constant AGGREGATION_QUERY_WITHOUT_ABNORMAL_HETEROZYGOSITY. */
	static final public String AGGREGATION_QUERY_WITHOUT_ABNORMAL_HETEROZYGOSITY = "WITHOUT_ABNORMAL_HETEROZYGOSITY";

	/** The Constant GENOTYPE_CODE_LABEL_ALL. */
	static final public String GENOTYPE_CODE_LABEL_ALL = "Any";

	/** The Constant GENOTYPE_CODE_LABEL_NOT_ALL_SAME. */
	static final public String GENOTYPE_CODE_LABEL_NOT_ALL_SAME = "Not all the same";

	/** The Constant GENOTYPE_CODE_LABEL_MOSTLY_SAME. */
	static final public String GENOTYPE_CODE_LABEL_MOSTLY_SAME = "All or mostly the same";

	/** The Constant GENOTYPE_CODE_LABEL_ALL_DIFFERENT. */
	static final public String GENOTYPE_CODE_LABEL_ALL_DIFFERENT = "All different";

	/** The Constant GENOTYPE_CODE_LABEL_NOT_ALL_DIFFERENT. */
	static final public String GENOTYPE_CODE_LABEL_NOT_ALL_DIFFERENT = "Not all different";

	/** The Constant GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_REF. */
	static final public String GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_REF = "All Homozygous Ref";

	/** The Constant GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_REF. */
	static final public String GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_REF = "Some Homozygous Ref";

	/** The Constant GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_VAR. */
	static final public String GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_VAR = "All Homozygous Var";

	/** The Constant GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_VAR. */
	static final public String GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_VAR = "Some Homozygous Var";

	/** The Constant GENOTYPE_CODE_LABEL_ALL_HETEROZYGOUS. */
	static final public String GENOTYPE_CODE_LABEL_ALL_HETEROZYGOUS = "All Heterozygous";

	/** The Constant GENOTYPE_CODE_LABEL_ATL_ONE_HETEROZYGOUS. */
	static final public String GENOTYPE_CODE_LABEL_ATL_ONE_HETEROZYGOUS = "Some Heterozygous";

	/** The Constant GENOTYPE_CODE_LABEL_WITHOUT_ABNORMAL_HETEROZYGOSITY. */
	static final public String GENOTYPE_CODE_LABEL_WITHOUT_ABNORMAL_HETEROZYGOSITY = "Without abnormal heterozygosity";

	/** The Constant genotypePatternToDescriptionMap. */
	static final private HashMap<String, String> genotypePatternToDescriptionMap = new LinkedHashMap<String, String>();

	/** The Constant genotypePatternToQueryMap. */
	static final private HashMap<String, String> genotypePatternToQueryMap = new HashMap<String, String>();

	static public final String MAIN_RESULT_PROJECTION_FIELD = "r";
	static public final String MOSTLY_SAME_GENOTYPE_RESULT_PROJECTION_FIELD = "r2";

	static
	{
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_ALL, "This will return all variants whithout applying any filters");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_NOT_ALL_SAME, "This will return variants where not all selected individuals have the same genotype");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_MOSTLY_SAME, "This will return variants where all or most selected individuals have the same genotype");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_ALL_DIFFERENT, "This will return variants where none of the selected individuals have the same genotype");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_NOT_ALL_DIFFERENT, "This will return variants where some of the selected individuals have the same genotypes");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_REF, "This will return variants where selected individuals are all homozygous with the reference allele");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_REF, "This will return variants where where at least one selected individual is homozygous with the reference allele");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_VAR, "This will return variants where selected individuals are all homozygous with an alternate allele");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_VAR, "This will return variants where at least one selected individual is homozygous with an alternate allele");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_ALL_HETEROZYGOUS, "This will return variants where selected individuals are all heterozygous");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_ATL_ONE_HETEROZYGOUS, "This will return variants where at least one selected individual is heterozygous");
		genotypePatternToDescriptionMap.put(GENOTYPE_CODE_LABEL_WITHOUT_ABNORMAL_HETEROZYGOSITY, "This will exclude variants where some alleles are found only in heterozygous genotypes (only for diploid, bi-allelic data)");
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_ALL, null);
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_MOSTLY_SAME, "$eq");
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_NOT_ALL_SAME, "$eq" + GenotypingDataQueryBuilder.AGGREGATION_QUERY_NEGATION_SUFFIX);
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_ALL_DIFFERENT, "$ne");
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_NOT_ALL_DIFFERENT, "$ne" + GenotypingDataQueryBuilder.AGGREGATION_QUERY_NEGATION_SUFFIX);
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_REF, "^0(/0)*$"/*|^$"*/ + GenotypingDataQueryBuilder.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX);
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_REF, "^0(/0)*$"/*|^$"*/ + GenotypingDataQueryBuilder.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX);
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_VAR, "^([1-9][0-9]*)(/\\1)*$"/*|^$"*/ + GenotypingDataQueryBuilder.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX);
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_VAR, "^([1-9][0-9]*)(/\\1)*$"/*|^$"*/ + GenotypingDataQueryBuilder.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX);
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_ALL_HETEROZYGOUS, "([0-9])([0-9])*(/(?!\\1))+([0-9])*"/*|^$"*/ + GenotypingDataQueryBuilder.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX);
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_ATL_ONE_HETEROZYGOUS, "([0-9])([0-9])*(/(?!\\1))+([0-9])*"/*|^$"*/ + GenotypingDataQueryBuilder.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX);
		genotypePatternToQueryMap.put(GENOTYPE_CODE_LABEL_WITHOUT_ABNORMAL_HETEROZYGOSITY, GenotypingDataQueryBuilder.AGGREGATION_QUERY_WITHOUT_ABNORMAL_HETEROZYGOSITY);
	}
	
	public GenotypingDataQueryBuilder(GigwaSearchVariantsRequest gsvr, MongoCollection<Document> tempExportColl, BasicDBList variantQueryDBList, boolean fForCounting) throws Exception
	{
		this.variantQueryDBList = variantQueryDBList;
		
        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

		this.mongoTemplate = MongoTemplateManager.get(sModule);

		this.genotypingProject = mongoTemplate.findById(Integer.valueOf(projId), GenotypingProject.class);
		this.geneNames = gsvr.getGeneName();
		this.variantEffects = gsvr.getVariantEffect();

		Query q = new Query();
       	q.addCriteria(Criteria.where("_id").is(projId));
       	q.addCriteria(Criteria.where(GenotypingProject.FIELDNAME_EFFECT_ANNOTATIONS + ".0").exists(true));
       	this.projectHasEffectAnnotations = mongoTemplate.findOne(q, GenotypingProject.class) != null;

       	this.selectedIndividuals[0] = gsvr.getCallSetIds().size() == 0 ? MgdbDao.getProjectIndividuals(sModule, projId) : gsvr.getCallSetIds().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(GigwaMethods.ID_SEPARATOR))).collect(Collectors.toSet());
		this.operator[0] = genotypePatternToQueryMap.get(gsvr.getGtPattern());
		this.mostSameRatio[0] = gsvr.getMostSameRatio();
		this.annotationFieldThresholds[0] = gsvr.getAnnotationFieldThresholds();
		this.missingData[0] = gsvr.getMissingData();
		this.minmaf[0] = gsvr.getMinmaf();
		this.maxmaf[0] = gsvr.getMaxmaf();
		final AtomicInteger nSampleCount = new AtomicInteger(0);
		this.individualToSampleListMap[0] = MgdbDao.getSamplesByIndividualForProject(sModule, projId, selectedIndividuals[0]);
		this.individualToSampleListMap[0].values().stream().map(spList -> nSampleCount.addAndGet(spList.size()));
		
		filteredGroups = getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gsvr, false);
		LOG.debug("Filtering genotypes on " + filteredGroups.size() + " groups");
		if (filteredGroups.contains(1))
		{
			this.selectedIndividuals[1] = gsvr.getCallSetIds2().size() == 0 ? MgdbDao.getProjectIndividuals(sModule, projId) : gsvr.getCallSetIds2().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(GigwaMethods.ID_SEPARATOR))).collect(Collectors.toSet());
			this.operator[1] = genotypePatternToQueryMap.get(gsvr.getGtPattern2());
			this.mostSameRatio[1] = gsvr.getMostSameRatio2();
			this.annotationFieldThresholds[1] = gsvr.getAnnotationFieldThresholds2();
			this.missingData[1] = gsvr.getMissingData2();
			this.minmaf[1] = gsvr.getMinmaf2();
			this.maxmaf[1] = gsvr.getMaxmaf2();
			this.individualToSampleListMap[1] = MgdbDao.getSamplesByIndividualForProject(sModule, projId, selectedIndividuals[1]);
			this.individualToSampleListMap[1].values().stream().map(spList -> nSampleCount.addAndGet(spList.size()));
			fDiscriminate = gsvr.isDiscriminate();
		}
		
		this.nTotalVariantCount = Helper.estimDocCount(mongoTemplate, MgdbDao.COLLECTION_NAME_TAGGED_VARIANT_IDS) + 1;
		if (this.nTotalVariantCount == 1)
		{
			MgdbDao.prepareDatabaseForSearches(mongoTemplate);	// list does not exist: create it
			this.nTotalVariantCount = Helper.estimDocCount(mongoTemplate,MgdbDao.COLLECTION_NAME_TAGGED_VARIANT_IDS) + 1;
		}
		this.taggedVariantList = mongoTemplate.findAll(Map.class, MgdbDao.COLLECTION_NAME_TAGGED_VARIANT_IDS);
		
		intervalIndexList = new ArrayList<>();
		for (int i=0; i<=taggedVariantList.size(); i++)
			intervalIndexList.add(i);
		
		this.groupFields = new Document();
		this.projectionFields = new Document();
		if (!fForCounting)
		{
	        groupFields.put(VariantData.FIELDNAME_REFERENCE_POSITION + "¤" + ReferencePosition.FIELDNAME_SEQUENCE, new Document("$first", "$" + VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE));
	        groupFields.put(VariantData.FIELDNAME_REFERENCE_POSITION + "¤" + ReferencePosition.FIELDNAME_START_SITE, new Document("$first", "$" + VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE));
	        groupFields.put(VariantData.FIELDNAME_REFERENCE_POSITION + "¤" + ReferencePosition.FIELDNAME_END_SITE, new Document("$first", "$" + VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_END_SITE));
	        groupFields.put(VariantData.FIELDNAME_TYPE, new Document("$first", "$" + VariantData.FIELDNAME_TYPE));
	        groupFields.put(VariantData.FIELDNAME_KNOWN_ALLELE_LIST, new Document("$first", "$" + VariantData.FIELDNAME_KNOWN_ALLELE_LIST));

	        boolean fIsMultiRunProject = genotypingProject.getRuns().size() > 1;
			projectionFields.put(AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, "$" + VariantData.FIELDNAME_REFERENCE_POSITION + (fIsMultiRunProject ? "¤" : ".") + ReferencePosition.FIELDNAME_SEQUENCE);
			projectionFields.put(AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, "$" + VariantData.FIELDNAME_REFERENCE_POSITION + (fIsMultiRunProject ? "¤" : ".") + ReferencePosition.FIELDNAME_START_SITE);
	        projectionFields.put(AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_END_SITE, "$" + VariantData.FIELDNAME_REFERENCE_POSITION + (fIsMultiRunProject ? "¤" : ".") + ReferencePosition.FIELDNAME_END_SITE);
	        projectionFields.put(AbstractVariantData.FIELDNAME_TYPE, "$" + VariantData.FIELDNAME_TYPE);
	        projectionFields.put(AbstractVariantData.FIELDNAME_KNOWN_ALLELE_LIST, "$" + VariantData.FIELDNAME_KNOWN_ALLELE_LIST);
	    }
    }
	
	/**
	 * Gets the number of queries.
	 *
	 * @return the number of queries
	 */
	public int getNumberOfQueries()
	{
		return (int) nTotalVariantCount;
	}

	/* (non-Javadoc)
	 * @see java.util.Iterator#remove()
	 */
	@Override
	public void remove()
	{
		throw new UnsupportedOperationException("Removal not supported");
	}
	
	/* (non-Javadoc)
	 * @see java.util.Iterator#hasNext()
	 */
	@Override
	public boolean hasNext()
	{
		return taggedVariantList.size() > 0 && intervalIndexList.size() > 0;
	}
	
	public List<Integer> suffleChunkOrder()
	{
		Collections.shuffle(intervalIndexList);
		return new ArrayList<Integer> (intervalIndexList);
	}

	/* (non-Javadoc)
	 * @see java.util.Iterator#next()
	 */
	@Override
	public List<BasicDBObject> next()
	{
		nNextCallCount++;
				        
		List<BasicDBObject> pipeline = new ArrayList<BasicDBObject>();
		BasicDBList initialMatchList = new BasicDBList(), annotationMatchList = new BasicDBList(), finalMatchList = new BasicDBList();

		pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", initialMatchList)));

		boolean fFilteringOnSequence = false;
		for (Object variantFilter : variantQueryDBList)
			if (((DBObject) variantFilter).containsField(AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE)) {
				fFilteringOnSequence = true;
				break;
			}
		
		if (fFilteringOnSequence)
			initialMatchList.addAll(variantQueryDBList);	// more efficient if added first in this case
		
        if (Helper.estimDocCount(mongoTemplate,GenotypingProject.class) != 1)
			initialMatchList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID, genotypingProject.getId()));
    	
        int currentInterval = intervalIndexList.get(0);
        intervalIndexList.remove(0);
        
    	BasicDBList chunkMatchAndList = new BasicDBList();
    	String leftBound = null, rightBound = null;
        if (currentInterval > 0) {
			leftBound = (String) taggedVariantList.get(currentInterval - 1).get("_id");
			chunkMatchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$gt", leftBound)));
		}
        
        if (currentInterval < taggedVariantList.size()) {
			rightBound = (String) taggedVariantList.get(currentInterval).get("_id");
			chunkMatchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$lte", rightBound)));
		}

		for (Object chunkMatch : chunkMatchAndList)
			initialMatchList.add((BasicDBObject) chunkMatch);

		/* Step to match variants according to annotations */			
		if (projectHasEffectAnnotations && (geneNames.length() > 0 || variantEffects.length() > 0)) {
			if (geneNames.length() > 0) {
				BasicDBObject geneNameDBO;
				if ("-".equals(geneNames))
					geneNameDBO = new BasicDBObject("$in", new String[] {"", null});
				else if ("+".equals(geneNames))
					geneNameDBO = new BasicDBObject("$regex", "^(?!\\s*$).+");
				else
					geneNameDBO = new BasicDBObject("$in", Helper.split(geneNames, ","));
				annotationMatchList.add(new BasicDBObject(VariantRunData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_GENE, geneNameDBO));
			}
			if (variantEffects.length() > 0)
				annotationMatchList.add(new BasicDBObject(VariantRunData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME, new BasicDBObject("$in", Helper.split(variantEffects, ","))));
			
			// we add another $match rather than re-using the first one because the first one may be pulled out from this pipeline and used for pre-filtering on the main variant collection, and functional annotations are only present in the run BasicDBObjects 
			pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", annotationMatchList)));
        }
		
        boolean fIsMultiRunProject = genotypingProject.getRuns().size() > 1;
		boolean[] fZygosityRegex = new boolean[2];
		boolean[] fIsWithoutAbnormalHeterozygosityQuery = new boolean[2];
		boolean[] fNegateMatch = new boolean[2];
		boolean[] fOr = new boolean[2];
		boolean[] fMafApplied = new boolean[2];
		boolean[] fMissingDataApplied = new boolean[2];
        boolean[] fCompareBetweenGenotypes = new boolean[2];
        String[] cleanOperator = new String[2];

        for (int g : filteredGroups) {
    		cleanOperator[g] = operator[g];
    		if (cleanOperator[g] != null) {
    			if (cleanOperator[g].endsWith(AGGREGATION_QUERY_NEGATION_SUFFIX)) {
    				fNegateMatch[g] = true;
    				cleanOperator[g] = cleanOperator[g].substring(0, cleanOperator[g].length() - AGGREGATION_QUERY_NEGATION_SUFFIX.length());
    			}
    			else if (cleanOperator[g].endsWith(AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX)) {
    				fZygosityRegex[g] = true;
    				cleanOperator[g] = cleanOperator[g].substring(0, cleanOperator[g].length() - AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX.length());
    			}
    			else if (cleanOperator[g].endsWith(AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX)) {
    				fZygosityRegex[g] = true;
    				fOr[g] = true;
    				cleanOperator[g] = cleanOperator[g].substring(0, cleanOperator[g].length() - AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX.length());
    			}
    			else if (cleanOperator[g].equals(AGGREGATION_QUERY_WITHOUT_ABNORMAL_HETEROZYGOSITY)) {
    				fIsWithoutAbnormalHeterozygosityQuery[g] = true;
    				cleanOperator[g] = cleanOperator[g].substring(0, cleanOperator[g].length() - AGGREGATION_QUERY_WITHOUT_ABNORMAL_HETEROZYGOSITY.length());
    			}
    		}

	        int nMaxNumberOfAllelesForOneVariant = maxAlleleCount > 0 ? maxAlleleCount : genotypingProject.getAlleleCounts().last(), nPloidy = genotypingProject.getPloidyLevel();
	        int nNumberOfPossibleGenotypes = (int) (nMaxNumberOfAllelesForOneVariant + MathUtils.factorial(nMaxNumberOfAllelesForOneVariant)/(MathUtils.factorial(nPloidy)*MathUtils.factorial(nMaxNumberOfAllelesForOneVariant-nPloidy)));
	        double maxMissingGenotypeCount = selectedIndividuals[g].size() * missingData[g] / 100;
    		if ("$ne".equals(cleanOperator[g]) && !fNegateMatch[g]) {
		        if (selectedIndividuals[g].size() - maxMissingGenotypeCount > nNumberOfPossibleGenotypes) {
		        	initialMatchList.add(new BasicDBObject("_id", null));	// return no results
		        	if (nNextCallCount == 1)
			        	LOG.info("Aborting 'all different' filter (more called individuals than possible genotypes in group " + (g + 1) + ")");
		        	return pipeline;
		        }
            }

			fCompareBetweenGenotypes[g] = cleanOperator[g] != null && !fZygosityRegex[g] && !fIsWithoutAbnormalHeterozygosityQuery[g];
	        if ("$ne".equals(cleanOperator[g]) && fNegateMatch[g]) {
		        if (selectedIndividuals[g].size() - maxMissingGenotypeCount > nNumberOfPossibleGenotypes) {
		        	fCompareBetweenGenotypes[g] = false;	// we know applying this filter would not affect the query
		        	if (nNextCallCount == 1)
			        	LOG.info("Ignoring 'not all different' filter on group 1 (more called individuals than possible genotypes in group " + (g + 1) + ")");
		        }
	        }
	        
	        fMafApplied[g] = maxmaf[g] != null && maxmaf[g].floatValue() < 50F || minmaf[g] != null && minmaf[g].floatValue() > 0.0F;
	        fMissingDataApplied[g] = missingData[g] != null && missingData[g] < 100;
        }
        
        if (variantQueryDBList.size() > 0 && !fFilteringOnSequence)
        	initialMatchList.addAll(variantQueryDBList);	// more efficient if added after chunking bit in this case
        
        BasicDBObject group = new BasicDBObject(groupFields);
		group.put("_id", "$_id." + VariantRunDataId.FIELDNAME_VARIANT_ID); // group multi-run records by variant id
		BasicDBObject project = new BasicDBObject(projectionFields);
        BasicDBObject addFieldsVars = new BasicDBObject();	// used for handling "all or mostly the same" filter
        BasicDBObject addFieldsIn = new BasicDBObject();	// used for handling "all or mostly the same" filter
    	BasicDBObject vars = new BasicDBObject();
    	BasicDBObject in = new BasicDBObject();
        BasicDBObject subIn = new BasicDBObject();
       
        for (int g : filteredGroups) {
    		boolean fMostSameSelected = "$eq".equals(cleanOperator[g]) && !fNegateMatch[g];
            boolean fNeedGtArray = fMissingDataApplied[g] || fMafApplied[g] || fZygosityRegex[g] || fIsWithoutAbnormalHeterozygosityQuery[g] || fCompareBetweenGenotypes[g] || fDiscriminate; // the only case when it's not needed is when we're only filtering on gene name or effect

            List<Object> currentGroupGtArray = new ArrayList<>();

            Iterator<String> indIt = selectedIndividuals[g].iterator();
            while (indIt.hasNext()) {
            	String ind = indIt.next();
				BasicDBList individualSampleGenotypeList = new BasicDBList();
				BasicDBList conditionsWhereAnnotationFieldValueIsTooLow = new BasicDBList();
				List<GenotypingSample> individualSamples = individualToSampleListMap[g].get(ind);
				for (int k=0; k<individualSamples.size(); k++) {	// this loop is executed only once for single-run projects
					GenotypingSample individualSample = individualSamples.get(k);
					String pathToGT = individualSample.getId() + "." + SampleGenotype.FIELDNAME_GENOTYPECODE;
					Object fullPathToGT = "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES/* + (int) ((individualSample.getId() - 1) / 100)*/ + "." + pathToGT;
					group.put(pathToGT.replaceAll("\\.", "¤"), new BasicDBObject("$addToSet", fullPathToGT));
					individualSampleGenotypeList.add("$" + pathToGT.replaceAll("\\.", "¤"));
	        		
					if (annotationFieldThresholds[g] != null)
						for (String annotation : annotationFieldThresholds[g].keySet()) {
							Float threshold = annotationFieldThresholds[g].get(annotation);
							if (threshold == 0)
								continue;

							String pathToAnnotationField = individualSample.getId() + "." + SampleGenotype.SECTION_ADDITIONAL_INFO + "." + annotation;
							group.put(pathToAnnotationField.replaceAll("\\.", "¤"), new BasicDBObject("$addToSet", "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + pathToAnnotationField));
							
							BasicDBList qualTooLowList = new BasicDBList();
							qualTooLowList.add(fIsMultiRunProject ? new BasicDBObject("$arrayElemAt", new Object[] {"$" + pathToAnnotationField.replaceAll("\\.", "¤"), 0}) : ("$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + pathToAnnotationField));
							qualTooLowList.add(threshold);
		
							BasicDBObject qualTooLow = new BasicDBObject("$lt", qualTooLowList);
							conditionsWhereAnnotationFieldValueIsTooLow.add(qualTooLow);
						}
	
			        if (k > 0)
						continue;	// the remaining code in this loop must only be executed once
			        
			        if (conditionsWhereAnnotationFieldValueIsTooLow.size() > 0)
			        	fullPathToGT = new BasicDBObject("$cond", new Object[] {new BasicDBObject("$or", conditionsWhereAnnotationFieldValueIsTooLow), null, fullPathToGT});
	                
	                if (fNeedGtArray && !fIsMultiRunProject)
	                	currentGroupGtArray.add(fullPathToGT);
		    	}
				if (fNeedGtArray && fIsMultiRunProject) {	// we're in the case of a multi-run project
					BasicDBObject union = new BasicDBObject("input", new BasicDBObject("$setUnion", individualSampleGenotypeList));
					union.put("as", "gt");
					union.put("cond", new BasicDBObject("$ne", Arrays.asList("$$gt", null)));
					BasicDBObject filteredGenotypeUnion = new BasicDBObject("$filter", union);	// union of (non-missing) genotypes for a given multi-sample individual
	
					currentGroupGtArray.add(conditionsWhereAnnotationFieldValueIsTooLow.size() == 0 ? filteredGenotypeUnion : new BasicDBObject("$cond", new Object[] { new BasicDBObject("$and", conditionsWhereAnnotationFieldValueIsTooLow), new Object[0], filteredGenotypeUnion}));
				}
			}

            if (fMafApplied[g]) {	// number of alternate alleles in selected population
            	BasicDBObject inObj;
				if (fIsMultiRunProject) {
	            	BasicDBList condList = new BasicDBList();
	                condList.add(new BasicDBObject("$eq", new Object[] {new BasicDBObject("$size", "$$g"), 1})); // if we have several distinct genotypes for this individual then we treat it as missing data (no alt allele to take into account)
	                condList.add(new BasicDBObject("$add", Arrays.asList(1, new BasicDBObject("$cmp", Arrays.asList(new BasicDBObject("$arrayElemAt", Arrays.asList("$$g", 0)), genotypingProject.getPloidyLevel() == 1 ? "1" : "0/1")))));
	                condList.add(0);
					inObj = new BasicDBObject("$cond", condList);
				}
				else
					inObj = new BasicDBObject("$add", Arrays.asList(1, new BasicDBObject("$cmp", Arrays.asList("$$g", genotypingProject.getPloidyLevel() == 1 ? "1" : "0/1"))));
            	in.put("a" + g, new BasicDBObject("$sum", new BasicDBObject("$map", new BasicDBObject("input", "$$gt" + g).append("as", "g").append("in", inObj))));
            }

            if (fMissingDataApplied[g] || fMafApplied[g] || (cleanOperator[g] != null && !fZygosityRegex[g]) || fIsWithoutAbnormalHeterozygosityQuery[g] || fDiscriminate) { //  number of missing genotypes in selected population
				if (fIsMultiRunProject)
					in.put("m" + g, new BasicDBObject("$sum", new BasicDBObject("$map", new BasicDBObject("input", "$$gt" + g).append("as", "g").append("in", new BasicDBObject("$abs", new BasicDBObject("$cmp", new Object[] {new BasicDBObject("$size", "$$g"), 1}))))));
				else// if (existingGenotypeCountList.size() > 0)
					in.put("m" + g, new BasicDBObject("$subtract", Arrays.asList(selectedIndividuals[g].size(), new BasicDBObject("$sum", new BasicDBObject("$map", new BasicDBObject("input", "$$gt" + g).append("as", "g").append("in", new BasicDBObject("$max", Arrays.asList(0, new BasicDBObject("$cmp", Arrays.asList("$$g", null))))))))));
            }
			
			if (fCompareBetweenGenotypes[g] && !fMostSameSelected) {
				BasicDBObject filter = new BasicDBObject("input", new BasicDBObject("$setUnion", "$$gt" + g));
				filter.put("as", "gt");
				filter.put("cond", fIsMultiRunProject ? new BasicDBObject("$eq", Arrays.asList(new BasicDBObject("$size", "$$gt"), 1)) : new BasicDBObject("$ne", Arrays.asList("$$gt", null)));
				in.put("dc" + g, new BasicDBObject("$size", new BasicDBObject("$filter", filter)));
			}
			else if (fZygosityRegex[g] || fIsWithoutAbnormalHeterozygosityQuery[g] || fMostSameSelected || fDiscriminate) {	//  distinct non-missing genotypes in selected population (zygosity comparison)
				if (fMostSameSelected || fDiscriminate)
					in.put("gt" + g, "$$gt" + g);	//  complete list of genotypes in selected population (all same)

				BasicDBObject filter = new BasicDBObject("input", new BasicDBObject("$setUnion", !fIsMultiRunProject ? "$$gt" + g : new BasicDBObject("$map", new BasicDBObject("input", "$$gt" + g).append("as", "g").append("in", new BasicDBObject("$arrayElemAt", Arrays.asList("$$g", 0))))));
				filter.put("as", "gt");
				filter.put("cond", new BasicDBObject("$ne", Arrays.asList("$$gt", null)));
				in.put("d" + g, new BasicDBObject("$filter", filter));
			}
			
			if (fMissingDataApplied[g])
				finalMatchList.add(new BasicDBObject("r.m" + g, new BasicDBObject("$lte", selectedIndividuals[g].size() * missingData[g] / 100)));
				
			if (fZygosityRegex[g] || fMissingDataApplied[g] || fMafApplied[g] || fCompareBetweenGenotypes[g] || fIsWithoutAbnormalHeterozygosityQuery[g] || fMostSameSelected || fDiscriminate) {	// we need to calculate extra fields via an additional $let operator
	            // keep previously computed fields
	            if (fMissingDataApplied[g] || (fCompareBetweenGenotypes[g]) || fDiscriminate)
	            	subIn.put("m" + g, "$$m" + g);
	            if (fZygosityRegex[g] || fMostSameSelected || fDiscriminate)
	            	subIn.put("d" + g, "$$d" + g);
	            if (fCompareBetweenGenotypes[g] && !fMostSameSelected)
	            	subIn.put("dc" + g, "$$dc" + g);

	            if (fCompareBetweenGenotypes[g] && !fMostSameSelected)
	            {	// dm = d + m
	            	 subIn.put("dm" + g, new BasicDBObject("$add", new Object[] {"$$dc" + g, "$$m" + g}));
	            	 
            		 finalMatchList.add(new BasicDBObject("r.m" + g, new BasicDBObject("$lt", selectedIndividuals[g].size() - 1)));	// if only one individual's genotype is not treated as missing then the filter makes no more sense
	            	 if ("$eq".equals(cleanOperator[g]) && fNegateMatch[g])
	           			 finalMatchList.add(new BasicDBObject("r.dc" + g, new BasicDBObject("$ne" /*not all same*/, 1)));
	            	 else if ("$ne".equals(cleanOperator[g]))
	            		 finalMatchList.add(new BasicDBObject("r.dm" + g, new BasicDBObject(fNegateMatch[g] ? "$lt" /*not all different*/ : "$eq" /*all different*/, selectedIndividuals[g].size())));
	            	 else
	            		 LOG.error("Invalid operator: " + operator);
	            }
	            
	            if (fMafApplied[g]) {	// allele frequency
	            	BasicDBList condList = new BasicDBList(), divideList = new BasicDBList();
	                condList.add(new BasicDBObject("$eq", new Object[] {"$$m" + g, selectedIndividuals[g].size()}));
	                condList.add(null);
	                condList.add(new BasicDBObject("$subtract", new Object[] {selectedIndividuals[g].size(), "$$m" + g}));
	                divideList.add(new BasicDBObject("$multiply", new Object[] {"$$a" + g, 50}));
	                divideList.add(new BasicDBObject("$cond", condList));

	        		subIn.put("f" + g, new BasicDBObject("$divide", divideList));
	        		
					BasicDBList orMafMatch = new BasicDBList();
					BasicDBList andMafMatch = new BasicDBList();
					andMafMatch.add(new BasicDBObject("r.f" + g, new BasicDBObject("$gte", minmaf[g])));
					andMafMatch.add(new BasicDBObject("r.f" + g, new BasicDBObject("$lte", maxmaf[g])));
					orMafMatch.add(new BasicDBObject("$and", andMafMatch));
					andMafMatch = new BasicDBList();
					andMafMatch.add(new BasicDBObject("r.f" + g, new BasicDBObject("$lte", Float.valueOf(100F - minmaf[g].floatValue()))));
					andMafMatch.add(new BasicDBObject("r.f" + g, new BasicDBObject("$gte", Float.valueOf(100F - maxmaf[g].floatValue()))));
					orMafMatch.add(new BasicDBObject("$and", andMafMatch));
					finalMatchList.add(new BasicDBObject("$or", orMafMatch));
	            }
	            
	            if (fIsWithoutAbnormalHeterozygosityQuery[g]) {	// counts for HZ, HR and HV genotypes
					BasicDBObject filter = new BasicDBObject("input", "$$d" + g);
					filter.put("as", "gt");
					filter.put("cond", new BasicDBObject("$eq", Arrays.asList("$$gt", "0/1")));
	            	subIn.put("hz" + g, new BasicDBObject("$size", new BasicDBObject("$filter", filter)));
	            	
	            	filter = new BasicDBObject("input", "$$d" + g);
					filter.put("as", "gt");
					filter.put("cond", new BasicDBObject("$eq", Arrays.asList("$$gt", "0/0")));
	            	subIn.put("hr" + g, new BasicDBObject("$size", new BasicDBObject("$filter", filter)));
					
	            	filter = new BasicDBObject("input", "$$d" + g);
					filter.put("as", "gt");
					filter.put("cond", new BasicDBObject("$eq", Arrays.asList("$$gt", "1/1")));
	            	subIn.put("hv" + g, new BasicDBObject("$size", new BasicDBObject("$filter", filter)));
	            }
	        }

			if (cleanOperator[g] != null || fDiscriminate) {
	            if (selectedIndividuals[g].size() >= 1) {
					if (fZygosityRegex[g]) {	// query to match specific genotype code with zygosity regex (homozygous var, homozygous ref, heterozygous)
						BasicDBList orSelectedGenotypeRegexAndFieldExistList = new BasicDBList();
						BasicDBObject orFinalSelectedGenotypeRegexAndFieldExist = new BasicDBObject();
						BasicDBObject andFinalSelectedGenotypeRegexAndFieldExist = new BasicDBObject();
					
						for (int j=0; j<selectedIndividuals[g].size(); j++) {
							/* FIXME: we can probably support heterozygous multiple-digit-genotypes using {$not : /^([0-9]+)(\/\1)*$/} */
							if (fOr[g]) {	// at least one whatever
								orSelectedGenotypeRegexAndFieldExistList.add(new BasicDBObject("r.d" + g + "." + j, new BasicDBObject("$regex", cleanOperator[g])));
							}
							else if (!genotypePatternToQueryMap.get(GENOTYPE_CODE_LABEL_ATL_ONE_HETEROZYGOUS).startsWith(cleanOperator[g])) {	// all homozygous whatever
								orSelectedGenotypeRegexAndFieldExistList.add(new BasicDBObject("r.d" + g + "." + j, new BasicDBObject(j == 0 ? "$regex" : "$exists", j == 0 ? cleanOperator[g] : false)));
							}
							else  {	// all heterozygous
								if (j == 0)
									orSelectedGenotypeRegexAndFieldExistList.add(new BasicDBObject("r.d" + g + "." + j, new BasicDBObject("$regex", cleanOperator[g])));
								else
								{
									BasicDBList orList = new BasicDBList();
									BasicDBObject clause1 = new BasicDBObject("r.d" + g + "." + j, new BasicDBObject("$exists", false));
									BasicDBObject clause2 = new BasicDBObject("r.d" + g + "." + j, new BasicDBObject("$regex", cleanOperator[g]));
									orList.add(clause1);
									orList.add(clause2);
									orSelectedGenotypeRegexAndFieldExistList.add(new BasicDBObject("$or", orList));
								}
							}
						}

						if (fOr[g]) {		
							orFinalSelectedGenotypeRegexAndFieldExist.put("$or", orSelectedGenotypeRegexAndFieldExistList);
							finalMatchList.add(orFinalSelectedGenotypeRegexAndFieldExist);
					    }			
						else {
							andFinalSelectedGenotypeRegexAndFieldExist.put("$and", orSelectedGenotypeRegexAndFieldExistList);
							finalMatchList.add(andFinalSelectedGenotypeRegexAndFieldExist);
						}
					}
					else if (fIsWithoutAbnormalHeterozygosityQuery[g]) {	// only for bi-allelic, diploid data: query that requires every allele present in heterozygous genotypes to be also present in homozygous ones
	                    BasicDBList orList = new BasicDBList();
	                    orList.add(new BasicDBObject("r.hz" + g, new BasicDBObject("$eq", 0)));
	                    BasicDBList andList = new BasicDBList();
	                    andList.add(new BasicDBObject("r.hv" + g, new BasicDBObject("$gt", 0)));
	                    andList.add(new BasicDBObject("r.hr" + g, new BasicDBObject("$gt", 0)));
	                    orList.add(new BasicDBObject("$and", andList));
	                    finalMatchList.add(new BasicDBObject("$or", orList));
	                }
					else if (fMostSameSelected || fDiscriminate) {
	                	BasicDBObject filter = new BasicDBObject("input", "$$gt" + g);
	                	filter.put("as", "g");
	                	filter.put("cond", new BasicDBObject("$eq", Arrays.asList("$$g", fIsMultiRunProject ? Arrays.asList("$$d") : "$$d")));
		                subIn.put("c" + g, new BasicDBObject("$map", new BasicDBObject("input", "$$d" + g).append("as", "d").append("in", new BasicDBObject("$size", new BasicDBObject("$filter", filter)))));
		                
		                addFieldsVars.put("dgc" + g, new BasicDBObject("$max", "$r.c" + g));	// dominant genotype count
		                Object minimumDominantGenotypeCount = new BasicDBObject("$multiply", Arrays.asList(new BasicDBObject("$subtract", new Object[] {selectedIndividuals[g].size(), "$r.m" + g}), mostSameRatio[g] / 100f));
		                
		                if (fMostSameSelected)
		                	addFieldsIn.put("ed" + g, new BasicDBObject("$gte", Arrays.asList("$$dgc" + g, minimumDominantGenotypeCount)));	// flag telling whether or not we have enough dominant genotypes to reach the required ratio
		                
		                if (fDiscriminate && g == 1) {
		                	addFieldsIn.put("dd", new BasicDBObject("$and", Arrays.asList(	/* dd (different dominant) set to true if both groups have exactly one dominant genotype and each group's dominant genotype differs from the other's */
	                			new BasicDBObject("$eq", Arrays.asList(1, new BasicDBObject("$size", new BasicDBObject("$filter", new BasicDBObject("input", "$r.c" + 0).append("cond", new BasicDBObject("$eq", Arrays.asList("$$dgc" + 0, "$$this"))))))),
	                			new BasicDBObject("$eq", Arrays.asList(1, new BasicDBObject("$size", new BasicDBObject("$filter", new BasicDBObject("input", "$r.c" + g).append("cond", new BasicDBObject("$eq", Arrays.asList("$$dgc" + g, "$$this"))))))),
	                			new BasicDBObject("$ne", Arrays.asList(new BasicDBObject("$arrayElemAt", Arrays.asList("$r.d" + 0, new BasicDBObject("$indexOfArray", Arrays.asList("$r.c" + 0, "$$dgc" + 0)))), new BasicDBObject("$arrayElemAt", Arrays.asList("$r.d" + g, new BasicDBObject("$indexOfArray", Arrays.asList("$r.c" + g, "$$dgc" + g))))))
	                		)));

		                	finalMatchList.add(new BasicDBObject("r2.dd", true));
		                }

		                if (fMostSameSelected)
		                	finalMatchList.add(new BasicDBObject("r2.ed" + g, true));
		            }
	            }
	        }

			vars.put("gt" + g, currentGroupGtArray);
		}
        
        if (subIn.size() > 0) { // insert additional $let
	        BasicDBObject subVars = in;
			BasicDBObject subLet = new BasicDBObject("vars", subVars);
			subLet.put("in", subIn);
			in = new BasicDBObject("$let", subLet);
        }

        if (in.size() > 0) {
			BasicDBObject let = new BasicDBObject("vars", vars);
			let.put("in", in);
			project.put(MAIN_RESULT_PROJECTION_FIELD, new BasicDBObject("$let", let));
        }
		
		if (fIsMultiRunProject)
			pipeline.add(new BasicDBObject("$group", group));
		if (!fIsMultiRunProject)
			project.put("_id", "$_id" + (!fIsMultiRunProject ? "." + VariantRunDataId.FIELDNAME_VARIANT_ID : ""));
		pipeline.add(new BasicDBObject("$project", project));
		
		if (addFieldsIn.size() > 0) {
			BasicDBObject addFieldsLet = new BasicDBObject("vars", addFieldsVars);
			addFieldsLet.put("in", addFieldsIn);
			pipeline.add(new BasicDBObject("$addFields", new BasicDBObject(MOSTLY_SAME_GENOTYPE_RESULT_PROJECTION_FIELD, new BasicDBObject("$let", addFieldsLet))));
		}
				
		if (finalMatchList.size() > 0)
			 pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", finalMatchList)));

//        if (nNextCallCount == 1) {
//        	try { System.out.println(new ObjectMapper().writeValueAsString(pipeline)); }
//        	catch (Exception ignored) {}
//        }
        return pipeline;
    }

//	private List<Comparable> buildFullListFromRange(MongoTemplate mongoTemplate, Comparable leftBound, Comparable rightBound) {
//		if (leftBound == null)
//			leftBound = mongoTemplate.findOne(new Query().with(new Sort(Sort.Direction.ASC, "_id")), VariantData.class).getId();
//		if (rightBound == null)
//			rightBound = mongoTemplate.findOne(new Query().with(new Sort(Sort.Direction.DESC, "_id")), VariantData.class).getId();
//
//		ArrayList<Comparable> result = new ArrayList<>();
//		String leftAsString = leftBound.toString(), rightAsString = rightBound.toString();
//		if (ObjectId.isValid(leftAsString) && ObjectId.isValid(rightAsString) && leftAsString.substring(0, 18).equals(rightAsString.substring(0, 18)))
//		{
//			int nCurrentId = Integer.parseInt(leftAsString.substring(18, 24), 16);
//			while (nCurrentId <= Integer.parseInt(rightAsString.substring(18, 24), 16))
//				result.add(new ObjectId(leftAsString.substring(0, 18) + Integer.toHexString(nCurrentId++)));
//		}
//		else
//		{
//			Query q = new Query().with(new Sort(Sort.Direction.ASC, "_id"));
//			q.fields().include("_id");
//			// not finished implementing method (functionality non needed)
//		}
//		return result;
//	}

	public static boolean areObjectIDsConsecutive(ObjectId first, ObjectId second)
	{
//		if (first == null || second == null)
//			return false;

		String firstAsString = first.toHexString(), secondAsString = second.toHexString();
		if (!firstAsString.substring(0, 18).equals(secondAsString.substring(0, 18)))
			return false;
		
		return 1 + Integer.parseInt(firstAsString.substring(18, 24), 16) == Integer.parseInt(secondAsString.substring(18, 24), 16);
	}

//	public static Document tryAndShrinkIdList(String pathToVariantId, Collection<Comparable> idCollection, int nShrinkThreshold)
//	{
//		if (idCollection.size() >= 300000)
//			try
//			{
//		//		long b4 = System.currentTimeMillis();
//		//		SortedSet<Comparable> idSet = SortedSet.class.isAssignableFrom(idCollection.getClass()) ? (SortedSet<Comparable>) idCollection : new TreeSet<Comparable>(idCollection);
//		//		System.out.println("sorting took " + (System.currentTimeMillis() - b4));
//				BasicDBList orList = new BasicDBList();
//				ArrayList<ObjectId> inIdList = new ArrayList<>(), rangeIdList = new ArrayList<>();
//				
//				ObjectId previousId = null;
//				for (Comparable id : idCollection)
//				{
//					ObjectId currentId = (ObjectId) id;
//					if (previousId == null || areObjectIDsConsecutive(previousId, currentId))
//						rangeIdList.add(currentId);
//					else
//					{
//						if (rangeIdList.size() >= nShrinkThreshold)
//						{	// replace list with a range
//							BasicDBList chunkMatchAndList = new BasicDBList();
//							chunkMatchAndList.add(new Document(pathToVariantId, new Document("$gte", rangeIdList.get(0))));
//							chunkMatchAndList.add(new Document(pathToVariantId, new Document("$lte", rangeIdList.get(rangeIdList.size() - 1))));
//							orList.add(new Document("$and", chunkMatchAndList));
//						}
//						else
//							inIdList.addAll(rangeIdList);	// range is too small, keep the list
//		
//						rangeIdList.clear();
//						rangeIdList.add(currentId);
//					}
//					previousId = currentId;
//				}
//				inIdList.addAll(rangeIdList);
//		
//				if (inIdList.size() > 0 || orList.size() == 0)
//					orList.add(new Document(pathToVariantId, new Document("$in", inIdList)));
//		
//				return orList.size() > 1 ? new Document("$or", orList) : (Document) orList.iterator().next();
//			}
//			catch (ClassCastException cce)
//			{
//				if (!cce.getMessage().contains("ObjectId"))
//					throw cce;	// otherwise it simply means IDs are of a different type, in which case we can't shrink the collection
//			}
////		else
////		{
////			LOG.debug("Didn't shrink id collection (" + idCollection.size() + " records only)");
////		}
//		
//		return new Document(pathToVariantId, new Document("$in", idCollection));	// not shrinked
//	}
	
	public static HashMap<String, String> getGenotypePatternToQueryMap() {
		return genotypePatternToQueryMap;
	}
	
	public static HashMap<String, String> getGenotypePatternToDescriptionMap() {
		return genotypePatternToDescriptionMap;
	}

	public void setMaxAlleleCount(int maxAlleleCount) {
		this.maxAlleleCount = maxAlleleCount;
	}

	static public List<Integer> getGroupsForWhichToFilterOnGenotypingData(GigwaSearchVariantsRequest gsvr, boolean fConsiderFieldThresholds)
    {
        List<Integer> result = new ArrayList<>();
    	if (gsvr.isDiscriminate() || !gsvr.getGtPattern().equals(GENOTYPE_CODE_LABEL_ALL) || (fConsiderFieldThresholds && gsvr.getAnnotationFieldThresholds().size() >= 1) || gsvr.getMissingData() < 100 || gsvr.getMinmaf() > 0 || gsvr.getMaxmaf() < 50)
    		result.add(0);
    	if (gsvr.isDiscriminate() || !gsvr.getGtPattern2().equals(GENOTYPE_CODE_LABEL_ALL) || (fConsiderFieldThresholds && gsvr.getAnnotationFieldThresholds2().size() >= 1) || gsvr.getMissingData2() < 100 || gsvr.getMinmaf2() > 0 || gsvr.getMaxmaf2() < 50)
    		result.add(1);

    	return result;
    }
	
	static public List<Integer> getGroupsForWhichToFilterOnGenotypingOrAnnotationData(GigwaSearchVariantsRequest gsvr, boolean fConsiderFielThresholds)
    {
        List<Integer> result = getGroupsForWhichToFilterOnGenotypingData(gsvr, fConsiderFielThresholds);
    	
    	if (result.size() == 0 && (Helper.estimDocCount(GigwaSearchVariantsRequest.getInfoFromId(gsvr.getVariantSetId(), 2)[0], GenotypingProject.class) != 1 || gsvr.getGeneName().length() > 0 || gsvr.getVariantEffect().length() > 0))
    		result.add(0);	// needed at least for filtering on annotation data or distinguish records according to project id

    	return result;
    }
}