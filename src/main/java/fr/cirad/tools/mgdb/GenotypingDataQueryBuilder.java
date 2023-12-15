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
import java.util.stream.Collectors;

import org.apache.commons.math.util.MathUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.bson.types.ObjectId;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;

import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongo.subtypes.VariantRunDataId;
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
    
    /** The operator. */
    private List<String> operator;
    
    private GigwaSearchVariantsRequest req;
    
    /** The individual index to sample list map. */
    private List<TreeMap<String /*individual*/, ArrayList<GenotypingSample>>> individualToSampleListMap = new ArrayList<>();
    
    /** The total chunk count. */
    private long nTotalChunkCount = 0;
    
    private int nNextCallCount = 0;
    
    private int maxAlleleCount = 0;
    
    private List<Integer> filteredGroups;
    
    /** Used for shuffling queried chunks */
    private List<Integer> intervalIndexList;
    
    private List<Map> taggedVariantList;

    private BasicDBList variantQueryDBList;
    
    boolean m_fFilteringOnSequence = false;

    private Document groupFields;

    private boolean fIsMultiRunProject = false;
    
    private boolean fForCounting = false;
    
    private boolean fGotMultiSampleIndividuals = false;
    
    private List<String> runsToRestrictQueryTo = null;
    
    private boolean fExcludeVariantsWithOnlyMissingData = false;
    
    private Document projectionFields = new Document();

    /** The Constant genotypePatternToDescriptionMap. */
    static final private HashMap<String, String> genotypePatternToDescriptionMap = new LinkedHashMap<String, String>();

    /** The Constant genotypePatternToQueryMap. */
    static final private HashMap<String, String> genotypePatternToQueryMap = new HashMap<String, String>();

    static public final String MAIN_RESULT_PROJECTION_FIELD = "r";

    static
    {
        genotypePatternToDescriptionMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL, "This will return all variants whithout applying any filters");
        genotypePatternToDescriptionMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_NOT_ALL_SAME, "This will return variants where not all selected individuals have the same genotype");
        genotypePatternToDescriptionMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_MOSTLY_SAME, "This will return variants where all or most selected individuals have the same genotype");
        genotypePatternToDescriptionMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_DIFFERENT, "This will return variants where none of the selected individuals have the same genotype");
        genotypePatternToDescriptionMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_NOT_ALL_DIFFERENT, "This will return variants where some of the selected individuals have the same genotypes");
        genotypePatternToDescriptionMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_REF, "This will return variants where selected individuals are all homozygous with the reference allele");
        genotypePatternToDescriptionMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_REF, "This will return variants where where at least one selected individual is homozygous with the reference allele");
        genotypePatternToDescriptionMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_VAR, "This will return variants where selected individuals are all homozygous with an alternate allele");
        genotypePatternToDescriptionMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_VAR, "This will return variants where at least one selected individual is homozygous with an alternate allele");
        genotypePatternToQueryMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL, null);
        genotypePatternToQueryMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_MOSTLY_SAME, "$eq");
        genotypePatternToQueryMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_NOT_ALL_SAME, "$eq" + GigwaSearchVariantsRequest.AGGREGATION_QUERY_NEGATION_SUFFIX);
        genotypePatternToQueryMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_DIFFERENT, "$ne");
        genotypePatternToQueryMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_NOT_ALL_DIFFERENT, "$ne" + GigwaSearchVariantsRequest.AGGREGATION_QUERY_NEGATION_SUFFIX);
        genotypePatternToQueryMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_REF, "^0(/0)*$"/*|^$"*/ + GigwaSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX);
        genotypePatternToQueryMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_REF, "^0(/0)*$"/*|^$"*/ + GigwaSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX);
        genotypePatternToQueryMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_VAR, "^([1-9][0-9]*)(/\\1)*$"/*|^$"*/ + GigwaSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX);
        genotypePatternToQueryMap.put(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_VAR, "^([1-9][0-9]*)(/\\1)*$"/*|^$"*/ + GigwaSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX);
    }

    public GenotypingDataQueryBuilder(GigwaSearchVariantsRequest gsvr, BasicDBList variantQueryDBList, boolean fForCounting) throws Exception
    {
        this.req = gsvr;
        this.fForCounting = fForCounting;
        this.variantQueryDBList = variantQueryDBList;
        Helper.convertIdFiltersToRunFormat(Arrays.asList(this.variantQueryDBList));
        
        for (Object variantFilter : variantQueryDBList)
            if (((BasicDBObject) variantFilter).containsKey(Assembly.getThreadBoundVariantRefPosPath() + "." + ReferencePosition.FIELDNAME_SEQUENCE)) {
                m_fFilteringOnSequence = true;
                break;
            }
        
        String info[] = Helper.getInfoFromId(req.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

        this.mongoTemplate = MongoTemplateManager.get(sModule);
        this.genotypingProject = mongoTemplate.findById(Integer.valueOf(projId), GenotypingProject.class);

        Query q = new Query();
        q.addCriteria(Criteria.where("_id").is(projId));
        q.addCriteria(Criteria.where(GenotypingProject.FIELDNAME_EFFECT_ANNOTATIONS + ".0").exists(true));
        this.projectHasEffectAnnotations = mongoTemplate.findOne(q, GenotypingProject.class) != null;

        filteredGroups = getGroupsForWhichToFilterOnGenotypingOrAnnotationData(req, false);
        int numberGroups = req.getNumberGroups();
        this.operator = new ArrayList<>(Collections.nCopies(numberGroups, null));
        this.individualToSampleListMap = new ArrayList<>(Collections.nCopies(numberGroups, new TreeMap<>()));

        LOG.debug("Filtering genotypes on " + filteredGroups.size() + " groups");
        List<List<String>> callsetIds = req.getAllCallSetIds();
        for (int nGroupIndex = 0; nGroupIndex < callsetIds.size(); nGroupIndex++) {
            Collection<String> groupIndividuals = callsetIds.isEmpty() || callsetIds.get(nGroupIndex).isEmpty() ? MgdbDao.getProjectIndividuals(sModule, projId) : callsetIds.get(nGroupIndex).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet());
            this.operator.set(nGroupIndex, genotypePatternToQueryMap.get(req.getGtPattern(nGroupIndex)));
            this.individualToSampleListMap.set(nGroupIndex, MgdbDao.getSamplesByIndividualForProject(sModule, projId, groupIndividuals));
        }
        
        this.nTotalChunkCount = Helper.estimDocCount(mongoTemplate, MgdbDao.COLLECTION_NAME_TAGGED_VARIANT_IDS) + 1;
        if (this.nTotalChunkCount == 1)
        {
            MgdbDao.prepareDatabaseForSearches(sModule);    // list does not exist: create it
            this.nTotalChunkCount = Helper.estimDocCount(mongoTemplate,MgdbDao.COLLECTION_NAME_TAGGED_VARIANT_IDS) + 1;
        }
        this.taggedVariantList = mongoTemplate.findAll(Map.class, MgdbDao.COLLECTION_NAME_TAGGED_VARIANT_IDS);

        intervalIndexList = new ArrayList<>();
        for (int i=0; i<=taggedVariantList.size(); i++)
            intervalIndexList.add(i);

        fIsMultiRunProject = genotypingProject.getRuns().size() > 1;

        for (TreeMap<String, ArrayList<GenotypingSample>> groupIndSamples : individualToSampleListMap)
            if (groupIndSamples != null && groupIndSamples.values().stream().filter(spList -> spList.size() > 1).findFirst().isPresent())
                fGotMultiSampleIndividuals = true;

        if (!fForCounting || fIsMultiRunProject) {
        	if (!individualToSampleListMap.isEmpty()) {
	            List<GenotypingSample> involvedSamples = new ArrayList<>();
	            for (int filteredGroup : filteredGroups) 
	                involvedSamples.addAll(individualToSampleListMap.get(filteredGroup).values().stream().flatMap(List::stream).collect(Collectors.toList()));
	
	            List<String> involvedProjectRuns = Helper.getRunsByProjectInSampleCollection(involvedSamples).get(genotypingProject.getId());
	            if (involvedProjectRuns != null && involvedProjectRuns.size() < genotypingProject.getRuns().size()) {
	                runsToRestrictQueryTo = involvedProjectRuns; // not all project runs are involved: adding a filter on the run field will make queries faster
	                fExcludeVariantsWithOnlyMissingData = true;  // some variants may have no data for the selected samples, we don't want to include them
	            }
        	}
            
            Integer nAssemblyId = Assembly.getThreadBoundAssembly();
            String refPosField = nAssemblyId != null ? AbstractVariantData.FIELDNAME_POSITIONS : AbstractVariantData.FIELDNAME_REFERENCE_POSITION;
            if (fIsMultiRunProject) {
                groupFields = new Document();
                Document firstContents = new Document();
                firstContents.append(ReferencePosition.FIELDNAME_SEQUENCE, "$" + refPosField + (nAssemblyId != null ? "." + nAssemblyId : "") + "." + ReferencePosition.FIELDNAME_SEQUENCE)
                            .append(ReferencePosition.FIELDNAME_START_SITE, "$" + refPosField + (nAssemblyId != null ? "." + nAssemblyId : "") + "." + ReferencePosition.FIELDNAME_START_SITE)
                            .append(ReferencePosition.FIELDNAME_END_SITE, "$" + refPosField + (nAssemblyId != null ? "." + nAssemblyId : "") + "." + ReferencePosition.FIELDNAME_END_SITE);
                groupFields.put(refPosField, new Document("$first", nAssemblyId != null ? new Document(nAssemblyId.toString(), firstContents) : firstContents));
                groupFields.put(VariantData.FIELDNAME_TYPE, new Document("$first", "$" + VariantData.FIELDNAME_TYPE));
                groupFields.put(VariantData.FIELDNAME_KNOWN_ALLELES, new Document("$first", "$" + VariantData.FIELDNAME_KNOWN_ALLELES));
            }
        }
    }
    
    /**
     * Gets the number of queries.
     *
     * @return the number of queries
     */
    public int getNumberOfQueries()
    {
        return (int) nTotalChunkCount;
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
    
    public List<Integer> shuffleChunkOrder()
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
        if (m_fFilteringOnSequence)
            initialMatchList.addAll(variantQueryDBList);    // more efficient if added first in this case

        int currentInterval = intervalIndexList.get(0);
        intervalIndexList.remove(0);
        
        BasicDBObject idRangeFilter = new BasicDBObject();
        if (currentInterval > 0)
        	idRangeFilter.append("$gt", taggedVariantList.get(currentInterval - 1).get("_id"));
        if (currentInterval < taggedVariantList.size())
        	idRangeFilter.append("$lte", taggedVariantList.get(currentInterval).get("_id"));
        if (!idRangeFilter.isEmpty())
            initialMatchList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, idRangeFilter));

        /* Step to match variants according to annotations */            
        if (projectHasEffectAnnotations && (req.getGeneName().length() > 0 || req.getVariantEffect().length() > 0)) {
            if (req.getGeneName().length() > 0) {
                BasicDBObject geneNameDBO;
                if ("-".equals(req.getGeneName()))
                    geneNameDBO = new BasicDBObject("$in", new String[] {"", null});
                else if ("+".equals(req.getGeneName()))
                    geneNameDBO = new BasicDBObject("$regex", "^(?!\\s*$).+");
                else
                    geneNameDBO = new BasicDBObject("$in", Helper.split(req.getGeneName(), ","));
                annotationMatchList.add(new BasicDBObject(VariantRunData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_GENE, geneNameDBO));
            }
            if (req.getVariantEffect().length() > 0)
                annotationMatchList.add(new BasicDBObject(VariantRunData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME, new BasicDBObject("$in", Helper.split(req.getVariantEffect(), ","))));
        }

        if (runsToRestrictQueryTo != null) {
            BasicDBList orList = new BasicDBList();
            orList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_RUNNAME, new BasicDBObject("$in", runsToRestrictQueryTo)));

            if (!annotationMatchList.isEmpty())
                orList.add(new BasicDBObject("$and", annotationMatchList));  // we do the annotation match here first, in case they would be in different run records than those that have the genotypes we want (FIXME: functional annotations should go into VariantData)

            initialMatchList.add(orList.size() == 1 ? orList.get(0) : new BasicDBObject("$or", orList));
        }
        pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", initialMatchList)));

        int maxi = filteredGroups.stream().max(Integer::compare).get() + 1;
        boolean[] fZygosityRegex = new boolean[maxi];
        boolean[] fNegateMatch = new boolean[maxi];
        boolean[] fOr = new boolean[maxi];
        boolean[] fMafApplied = new boolean[maxi];
        boolean[] fMissingDataApplied = new boolean[maxi];
        boolean[] fHezRatioApplied = new boolean[maxi];
        boolean[] fCompareBetweenGenotypes = new boolean[maxi];
        String[] cleanOperator = new String[maxi];

        if (req.getNumberGroups() > 0)
            for (int g : filteredGroups) {
                cleanOperator[g] = operator.get(g);
                if (cleanOperator[g] != null) {
                    if (cleanOperator[g].endsWith(GigwaSearchVariantsRequest.AGGREGATION_QUERY_NEGATION_SUFFIX)) {
                        fNegateMatch[g] = true;
                        cleanOperator[g] = cleanOperator[g].substring(0, cleanOperator[g].length() - GigwaSearchVariantsRequest.AGGREGATION_QUERY_NEGATION_SUFFIX.length());
                    }
                    else if (cleanOperator[g].endsWith(GigwaSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX)) {
                        fZygosityRegex[g] = true;
                        cleanOperator[g] = cleanOperator[g].substring(0, cleanOperator[g].length() - GigwaSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX.length());
                    }
                    else if (cleanOperator[g].endsWith(GigwaSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX)) {
                        fZygosityRegex[g] = true;
                        fOr[g] = true;
                        cleanOperator[g] = cleanOperator[g].substring(0, cleanOperator[g].length() - GigwaSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX.length());
                    }
                }
    
                int nMaxNumberOfAllelesForOneVariant = maxAlleleCount > 0 ? maxAlleleCount : genotypingProject.getAlleleCounts().last(), nPloidy = genotypingProject.getPloidyLevel();
                int nNumberOfPossibleGenotypes;
                try {
                    nNumberOfPossibleGenotypes = (int) (nMaxNumberOfAllelesForOneVariant + MathUtils.factorial(nMaxNumberOfAllelesForOneVariant)/(MathUtils.factorial(nPloidy)*MathUtils.factorial(nMaxNumberOfAllelesForOneVariant-nPloidy)));
                }
                catch (ArithmeticException ae) {    // nMaxNumberOfAllelesForOneVariant must be too large for its factorial to be calculated
                    nNumberOfPossibleGenotypes = Integer.MAX_VALUE;
                }
                
                int nGroupSize = individualToSampleListMap.get(g).size();
                double maxMissingGenotypeCount = nGroupSize * req.getMaxMissingData(g) / 100;
                if ("$ne".equals(cleanOperator[g]) && !fNegateMatch[g]) {
                    if (nGroupSize - maxMissingGenotypeCount > nNumberOfPossibleGenotypes) {
                        initialMatchList.add(new BasicDBObject("_id", null));    // return no results
                        if (nNextCallCount == 1)
                            LOG.info("Aborting 'all different' filter (more called individuals than possible genotypes in group " + (g + 1) + ")");
                        return pipeline;
                    }
                }
    
                fCompareBetweenGenotypes[g] = cleanOperator[g] != null && !fZygosityRegex[g];
                if ("$ne".equals(cleanOperator[g]) && fNegateMatch[g]) {
                    if (nGroupSize - maxMissingGenotypeCount > nNumberOfPossibleGenotypes) {
                        fCompareBetweenGenotypes[g] = false;    // we know applying this filter would not affect the query
                        if (nNextCallCount == 1)
                            LOG.info("Ignoring 'not all different' filter on group 1 (more called individuals than possible genotypes in group " + (g + 1) + ")");
                    }
                }
                
                fMafApplied[g] = req.getMaxMaf(g) != null && req.getMaxMaf(g) < 50F || req.getMinMaf(g) != null && req.getMinMaf(g) > 0.0F;
                fMissingDataApplied[g] = (req.getMinMissingData(g) != null && req.getMinMissingData(g) > 0) || (req.getMaxMissingData(g) != null && req.getMaxMissingData(g) < 100);
                fHezRatioApplied[g] = (req.getMinHeZ(g) != null && req.getMinHeZ(g) > 0) || (req.getMaxHeZ(g) != null && req.getMaxHeZ(g) < 100);
                
                if (req.getMinHeZ(g) != null && req.getMinHeZ(g) > 0 && fZygosityRegex[g]) {
                    if (fOr[g]) {   // check if it's possible to have any homozygous genotypes based on selected heterozygosity and missing data rates
                        double minMissing = !fMissingDataApplied[g] ? 0 : Math.ceil(req.getMinMissingData(g) * nGroupSize / 100);
                        double minNonHomoZ = Math.ceil((req.getMinHeZ(g)) * (nGroupSize * (1 - (!fMissingDataApplied[g] ? 0 : req.getMinMissingData(g)) / 100f)) / 100);
                        if (minNonHomoZ + minMissing >= nGroupSize) {
                            initialMatchList.add(new BasicDBObject("_id", null));    // return no results
                            if (nNextCallCount == 1)
                                LOG.info("Aborting incompatible heterozygous / homozygous filters (no space for any homozygous)");
                            return pipeline;
                        }
                    }
                    else {
                        initialMatchList.add(new BasicDBObject("_id", null));    // return no results
                        if (nNextCallCount == 1)
                            LOG.info("Aborting incompatible heterozygous / homozygous filters (no space for any heterozygous)");
                        return pipeline;
                    }
                }
            }

        if (variantQueryDBList.size() > 0 && !m_fFilteringOnSequence)
            initialMatchList.addAll(variantQueryDBList);    // more efficient if added after chunking bit in this case
        
        if (fIsMultiRunProject)
            groupFields.put("_id", "$_id." + VariantRunDataId.FIELDNAME_VARIANT_ID); // group multi-run records by variant id

        BasicDBObject addFieldsVars = new BasicDBObject();    // used for handling "all or mostly the same" filter
        BasicDBObject addFieldsIn = new BasicDBObject();    // used for handling "all or mostly the same" filter
        BasicDBObject vars = new BasicDBObject();
        BasicDBObject in = new BasicDBObject();
        BasicDBObject subIn = new BasicDBObject();
       
        if (req.getNumberGroups() > 0)
            for (int g : filteredGroups) {
                boolean fMostSameSelected = "$eq".equals(cleanOperator[g]) && !fNegateMatch[g];
                boolean fNeedGtArray = fExcludeVariantsWithOnlyMissingData || fMissingDataApplied[g] || fMafApplied[g] || fZygosityRegex[g] || fHezRatioApplied[g] || fCompareBetweenGenotypes[g] || req.isDiscriminate(g); // the only case when it's not needed is when we're only filtering on gene name or effect
                int nGroupSize = individualToSampleListMap.get(g).size();
                
                List<Object> currentGroupGtArray = new ArrayList<>();
                for (String ind : individualToSampleListMap.get(g).keySet()) {
                    BasicDBList individualSampleGenotypeList = new BasicDBList();
                    BasicDBList conditionsWhereAnnotationFieldValueIsTooLow = new BasicDBList();
                    List<GenotypingSample> individualSamples = individualToSampleListMap.get(g).get(ind);
                    if (individualSamples != null)
                        for (int k=0; k<individualSamples.size(); k++) {    // this loop is executed only once for single-run projects
                            GenotypingSample individualSample = individualSamples.get(k);
                            Object fullPathToGT = "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + individualSample.getId() + "." + SampleGenotype.FIELDNAME_GENOTYPECODE;
                            if (fNeedGtArray && fIsMultiRunProject)
                                groupFields.put(SampleGenotype.FIELDNAME_GENOTYPECODE + "_" + individualSample.getId(), new BasicDBObject("$addToSet", fullPathToGT));
                            individualSampleGenotypeList.add("$" + SampleGenotype.FIELDNAME_GENOTYPECODE + "_" + individualSample.getId());
                            
                            if (req.getAnnotationFieldThresholds(g) != null)
                                for (String annotation : req.getAnnotationFieldThresholds(g).keySet()) {
                                    Float threshold = req.getAnnotationFieldThresholds(g).get(annotation);
                                    if (threshold == 0)
                                        continue;
        
                                    String pathToAnnotationField = individualSample.getId() + "." + SampleGenotype.SECTION_ADDITIONAL_INFO + "." + annotation;
                                    if (fNeedGtArray && fIsMultiRunProject)
                                        groupFields.put(pathToAnnotationField, new BasicDBObject("$addToSet", "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + pathToAnnotationField));
                                    
                                    BasicDBList qualTooLowList = new BasicDBList();
                                    qualTooLowList.add(fIsMultiRunProject ? new BasicDBObject("$arrayElemAt", new Object[] {"$" + pathToAnnotationField, 0}) : ("$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + pathToAnnotationField));
                                    qualTooLowList.add(threshold);
                
                                    BasicDBObject qualTooLow = new BasicDBObject("$lt", qualTooLowList);
                                    conditionsWhereAnnotationFieldValueIsTooLow.add(qualTooLow);
                                }
            
                            if (k > 0)
                                continue;    // the remaining code in this loop must only be executed once
                            
                            if (conditionsWhereAnnotationFieldValueIsTooLow.size() > 0)
                                fullPathToGT = new BasicDBObject("$cond", new Object[] {new BasicDBObject("$or", conditionsWhereAnnotationFieldValueIsTooLow), null, fullPathToGT});
                            
                            if (fNeedGtArray && !fIsMultiRunProject)
                                currentGroupGtArray.add(fullPathToGT);
                        }
                    if (fNeedGtArray && fIsMultiRunProject) {    // we're in the case of a multi-run project
                        BasicDBObject union = new BasicDBObject("input", new BasicDBObject("$setUnion", individualSampleGenotypeList));
                        union.put("as", "gt");
                        union.put("cond", new BasicDBObject("$ne", Arrays.asList("$$gt", null)));
                        BasicDBObject filteredGenotypeUnion = new BasicDBObject("$filter", union);    // union of (non-missing) genotypes for a given multi-sample individual
        
                        currentGroupGtArray.add(conditionsWhereAnnotationFieldValueIsTooLow.size() == 0 ? filteredGenotypeUnion : new BasicDBObject("$cond", new Object[] { new BasicDBObject("$and", conditionsWhereAnnotationFieldValueIsTooLow), new Object[0], filteredGenotypeUnion}));
                    }
                }
    
                if (fMafApplied[g]) {    // number of alternate alleles in selected population
                    BasicDBObject inObj;
                    if (fIsMultiRunProject) {
                        BasicDBList condList = new BasicDBList();
                        BasicDBObject addObject = new BasicDBObject("$add", Arrays.asList(1, new BasicDBObject("$cmp", Arrays.asList(new BasicDBObject("$arrayElemAt", Arrays.asList("$$g", 0)), genotypingProject.getPloidyLevel() == 1 ? "1" : "0/1"))));
                        if (!fGotMultiSampleIndividuals)
                            inObj = addObject;    // no need to make sure all genotypes for each individual are equal because there's only one sample per individual 
                        else {
                            condList.add(new BasicDBObject("$eq", new Object[] {new BasicDBObject("$size", "$$g"), 1})); // if we have several distinct genotypes for this individual then we treat it as missing data (no alt allele to take into account)
                            condList.add(addObject);
                            condList.add(0);
                            inObj = new BasicDBObject("$cond", condList);
                        }
                    }
                    else
                        inObj = new BasicDBObject("$add", Arrays.asList(1, new BasicDBObject("$cmp", Arrays.asList("$$g", genotypingProject.getPloidyLevel() == 1 ? "1" : "0/1"))));
                    in.put("a" + g, new BasicDBObject("$sum", new BasicDBObject("$map", new BasicDBObject("input", "$$gt" + g).append("as", "g").append("in", inObj))));
                }
    
                if (fExcludeVariantsWithOnlyMissingData || fMissingDataApplied[g] || fMafApplied[g] || fHezRatioApplied[g] || (cleanOperator[g] != null && !fZygosityRegex[g]) || req.isDiscriminate(g)) { //  number of missing genotypes in selected population
                    if (fIsMultiRunProject)
                        in.put("m" + g, new BasicDBObject("$sum", new BasicDBObject("$map", new BasicDBObject("input", "$$gt" + g).append("as", "g").append("in", new BasicDBObject("$abs", new BasicDBObject("$cmp", new Object[] {new BasicDBObject("$size", "$$g"), 1}))))));
                    else// if (existingGenotypeCountList.size() > 0)
                        in.put("m" + g, new BasicDBObject("$subtract", Arrays.asList(nGroupSize, new BasicDBObject("$sum", new BasicDBObject("$map", new BasicDBObject("input", "$$gt" + g).append("as", "g").append("in", new BasicDBObject("$max", Arrays.asList(0, new BasicDBObject("$cmp", Arrays.asList("$$g", null))))))))));
                }
    
                if (fCompareBetweenGenotypes[g] && !fMostSameSelected) {
                    BasicDBObject filter = new BasicDBObject("input", new BasicDBObject("$setUnion", "$$gt" + g));
                    filter.put("as", "gt");
                    filter.put("cond", fIsMultiRunProject ? new BasicDBObject("$eq", Arrays.asList(new BasicDBObject("$size", "$$gt"), 1)) : new BasicDBObject("$ne", Arrays.asList("$$gt", null)));
                    in.put("dc" + g, new BasicDBObject("$size", new BasicDBObject("$filter", filter)));
                }
    
                if (fZygosityRegex[g]|| fMostSameSelected || req.isDiscriminate(g)) {    //  distinct non-missing genotypes in selected population (zygosity comparison)
                    if (fMostSameSelected || req.isDiscriminate(g))
                        in.put("gt" + g, "$$gt" + g);    //  complete list of genotypes in selected population (all same)
    
                    BasicDBObject filter = new BasicDBObject("input", new BasicDBObject("$setUnion", !fIsMultiRunProject ? "$$gt" + g : new BasicDBObject("$map", new BasicDBObject("input", "$$gt" + g).append("as", "g").append("in", new BasicDBObject("$arrayElemAt", Arrays.asList("$$g", 0))))));
                    filter.put("as", "gt");
                    filter.put("cond", new BasicDBObject("$ne", Arrays.asList("$$gt", null)));
                    in.put("d" + g, new BasicDBObject("$filter", filter));
                }
    
                if (fMissingDataApplied[g]) {
                    BasicDBObject missingDataFilter = new BasicDBObject();
                    if (req.getMinMissingData(g) != null && req.getMinMissingData(g) > 0)
                        missingDataFilter.put("$gte", nGroupSize * req.getMinMissingData(g) / 100);
                    if (req.getMaxMissingData(g) != null && req.getMaxMissingData(g) < 100)
                        missingDataFilter.put("$lte", nGroupSize * req.getMaxMissingData(g) / 100);
                    finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".m" + g, missingDataFilter));
                }
    
                if (fHezRatioApplied[g]) {    // heterozygosity ratio
                    BasicDBObject filter = new BasicDBObject("input", !fIsMultiRunProject ? "$$gt" + g : new BasicDBObject("$map", new BasicDBObject("input", "$$gt" + g).append("as", "g").append("in", new BasicDBObject("$arrayElemAt", Arrays.asList("$$g", 0)))));
                    filter.put("as", "gt");
                    filter.put("cond", new BasicDBObject("$and", Arrays.asList(new BasicDBObject("$ne", Arrays.asList("$$gt", null)), new BasicDBObject("$not", new BasicDBObject("$regexMatch", new BasicDBObject("input", "$$gt").append("regex", "^([0-9]+)(\\/\\1)*$"))))));
                    in.put("he" + g, new BasicDBObject("$size", new BasicDBObject("$filter", filter))); // heterozygous genotype count
    
                    BasicDBObject hzFilter = new BasicDBObject();
                    if (req.getMinHeZ(g) != null && req.getMinHeZ(g) > 0)
                        hzFilter.put("$gte", req.getMinHeZ(g));
                    if (req.getMaxHeZ(g) != null && req.getMaxHeZ(g) < 100)
                        hzFilter.put("$lte", req.getMaxHeZ(g));
                    finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".hef" + g, hzFilter));
                }
    
                if (fZygosityRegex[g] || fExcludeVariantsWithOnlyMissingData || fMissingDataApplied[g] || fMafApplied[g] || fCompareBetweenGenotypes[g] || fHezRatioApplied[g] || fMostSameSelected || req.isDiscriminate(g)) {    // we need to calculate extra fields via an additional $let operator
                    // keep previously computed fields
                    if (fExcludeVariantsWithOnlyMissingData || fMissingDataApplied[g] || (fCompareBetweenGenotypes[g]) || fHezRatioApplied[g] || req.isDiscriminate(g))
                        subIn.put("m" + g, "$$m" + g);
                    if (fZygosityRegex[g] || fMostSameSelected || req.isDiscriminate(g))
                        subIn.put("d" + g, "$$d" + g);
                    if (fCompareBetweenGenotypes[g] && !fMostSameSelected)
                        subIn.put("dc" + g, "$$dc" + g);
                    if (fHezRatioApplied[g])
                        subIn.put("he" + g, "$$he" + g);
    
                    if (fCompareBetweenGenotypes[g] && !fMostSameSelected)
                    {    // dm = d + m
                         subIn.put("dm" + g, new BasicDBObject("$add", new Object[] {"$$dc" + g, "$$m" + g}));
                         
                         finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".m" + g, new BasicDBObject("$lt", nGroupSize - 1)));    // if only one individual's genotype is not treated as missing then the filter makes no more sense
                         if ("$eq".equals(cleanOperator[g]) && fNegateMatch[g])
                                finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".dc" + g, new BasicDBObject("$ne" /*not all same*/, 1)));
                         else if ("$ne".equals(cleanOperator[g]))
                             finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".dm" + g, new BasicDBObject(fNegateMatch[g] ? "$lt" /*not all different*/ : "$eq" /*all different*/, nGroupSize)));
                         else
                             LOG.error("Invalid operator: " + operator);
                    }
    
                    if (fMafApplied[g]) {    // allele frequency
                        BasicDBList condList = new BasicDBList(), divideList = new BasicDBList();
                        condList.add(new BasicDBObject("$eq", new Object[] {"$$m" + g, nGroupSize}));
                        condList.add(null);
                        condList.add(new BasicDBObject("$subtract", new Object[] {nGroupSize, "$$m" + g}));
                        divideList.add(new BasicDBObject("$multiply", new Object[] {"$$a" + g, 50}));
                        divideList.add(new BasicDBObject("$cond", condList));
    
                        subIn.put("f" + g, new BasicDBObject("$divide", divideList));
                        
                        float minMaf = req.getMinMaf(g), maxMaf = req.getMaxMaf(g);
                        BasicDBList orMafMatch = new BasicDBList();
                        BasicDBList andMafMatch = new BasicDBList();
                        if (minMaf > 0 && minMaf != maxMaf)
                            andMafMatch.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".f" + g, new BasicDBObject("$gte", minMaf)));
                        andMafMatch.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".f" + g, new BasicDBObject(minMaf == maxMaf ? "$eq" :"$lte", maxMaf != 50 ? maxMaf : (100F - minMaf))));
                        orMafMatch.add(andMafMatch.size() == 1 ? andMafMatch.iterator().next() : new BasicDBObject("$and", andMafMatch));
                        if (maxMaf != 50) {
                            andMafMatch = new BasicDBList();
                            if (minMaf > 0 && minMaf != maxMaf)
                                andMafMatch.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".f" + g, new BasicDBObject("$lte", (100F - minMaf))));
                            andMafMatch.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".f" + g, new BasicDBObject(minMaf == maxMaf ? "$eq" :"$gte", (100F - maxMaf))));
                            orMafMatch.add(andMafMatch.size() == 1 ? andMafMatch.iterator().next() : new BasicDBObject("$and", andMafMatch));
                        }
                        finalMatchList.add(new BasicDBObject("$or", orMafMatch));
                    }
                }
    
                if (cleanOperator[g] != null || req.isDiscriminate(g)) {
                    if (nGroupSize >= 1) {
                        if (fZygosityRegex[g]) {    // query to match specific genotype code with zygosity regex (homozygous var, homozygous ref, heterozygous)
                            BasicDBList atLeastOneSelectedGenotypeRegexAndFieldExistList = new BasicDBList();
                            BasicDBObject atLeastOneFinalSelectedGenotypeRegexAndFieldExist = new BasicDBObject();
                            BasicDBObject allFinalSelectedGenotypeRegexAndFieldExist = new BasicDBObject();
                        
                            for (int j = 0; j < nGroupSize; j++) {
                                if (fOr[g]) { // at least one homozygous
                                    atLeastOneSelectedGenotypeRegexAndFieldExistList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".d" + g + "." + j, new BasicDBObject("$regex", cleanOperator[g])));
                                }
                                else if (j <= 1) { // all homozygous
                                    atLeastOneSelectedGenotypeRegexAndFieldExistList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".d" + g + "." + j, new BasicDBObject(j == 0 ? "$regex" : "$exists", j == 0 ? cleanOperator[g] : false)));
                                }
                            }
    
                            if (fOr[g]) {        
                                atLeastOneFinalSelectedGenotypeRegexAndFieldExist.put("$or", atLeastOneSelectedGenotypeRegexAndFieldExistList);
                                finalMatchList.add(atLeastOneFinalSelectedGenotypeRegexAndFieldExist);
                            }            
                            else {
                                allFinalSelectedGenotypeRegexAndFieldExist.put("$and", atLeastOneSelectedGenotypeRegexAndFieldExistList);
                                finalMatchList.add(allFinalSelectedGenotypeRegexAndFieldExist);
                            }
                        }
    
                        if (fMostSameSelected || req.isDiscriminate(g)) {    
                            BasicDBObject filter = new BasicDBObject("input", "$$gt" + g);
                            filter.put("as", "g");
                            filter.put("cond", new BasicDBObject("$eq", Arrays.asList("$$g", fIsMultiRunProject ? Arrays.asList("$$d") : "$$d")));
                            subIn.put("c" + g, new BasicDBObject("$map", new BasicDBObject("input", "$$d" + g).append("as", "d").append("in", new BasicDBObject("$size", new BasicDBObject("$filter", filter)))));
                            
                            addFieldsVars.put("dgc" + g, new BasicDBObject("$max", "$" + MAIN_RESULT_PROJECTION_FIELD + ".c" + g));    // dominant genotype count
                            Object minimumDominantGenotypeCount = new BasicDBObject("$multiply", Arrays.asList(new BasicDBObject("$subtract", new Object[] {nGroupSize, "$" + MAIN_RESULT_PROJECTION_FIELD + ".m" + g}), req.getMostSameRatio(g) / 100f));
                            
                            if (fMostSameSelected)
                                addFieldsIn.put("ed" + g, new BasicDBObject("$gte", Arrays.asList("$$dgc" + g, minimumDominantGenotypeCount)));    // flag telling whether or not we have enough dominant genotypes to reach the required ratio
                            
                            if (fMissingDataApplied[g] || fMostSameSelected)
                                addFieldsIn.put("m" + g, "$" + MAIN_RESULT_PROJECTION_FIELD + ".m" + g);
                            
                            if (fMafApplied[g])
                                addFieldsIn.put("f" + g, "$" + MAIN_RESULT_PROJECTION_FIELD + ".f" + g); 
    
                            Integer groupToDiscriminateWith = req.getDiscriminate().get(g);
                            if (groupToDiscriminateWith != null) {
                            	groupToDiscriminateWith--;
                                addFieldsIn.put("dd" + g, new BasicDBObject("$and", Arrays.asList(    /* dd (different dominant) set to true if both groups have exactly one dominant genotype and each group's dominant genotype differs from the other's */
                                    new BasicDBObject("$eq", Arrays.asList(1, new BasicDBObject("$size", new BasicDBObject("$filter", new BasicDBObject("input", "$" + MAIN_RESULT_PROJECTION_FIELD + ".c" + groupToDiscriminateWith).append("cond", new BasicDBObject("$eq", Arrays.asList("$$dgc" + groupToDiscriminateWith, "$$this"))))))),
                                    new BasicDBObject("$eq", Arrays.asList(1, new BasicDBObject("$size", new BasicDBObject("$filter", new BasicDBObject("input", "$" + MAIN_RESULT_PROJECTION_FIELD + ".c" + g).append("cond", new BasicDBObject("$eq", Arrays.asList("$$dgc" + g, "$$this"))))))),
                                    new BasicDBObject("$ne", Arrays.asList(new BasicDBObject("$arrayElemAt", Arrays.asList("$" + MAIN_RESULT_PROJECTION_FIELD + ".d" + groupToDiscriminateWith, new BasicDBObject("$indexOfArray", Arrays.asList("$" + MAIN_RESULT_PROJECTION_FIELD + ".c" + groupToDiscriminateWith, "$$dgc" + groupToDiscriminateWith)))), new BasicDBObject("$arrayElemAt", Arrays.asList("$" + MAIN_RESULT_PROJECTION_FIELD + ".d" + g, new BasicDBObject("$indexOfArray", Arrays.asList("$" + MAIN_RESULT_PROJECTION_FIELD + ".c" + g, "$$dgc" + g))))))
                                )));
    
                                finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".dd" + g, true));
                            }
    
                            if (fMostSameSelected)
                                finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".ed" + g, true));
                        }
                    }
                }
                
                if (fHezRatioApplied[g]) {
                    BasicDBList condList = new BasicDBList();
                    condList.add(new BasicDBObject("$eq", new Object[] {"$" + MAIN_RESULT_PROJECTION_FIELD + ".m" + g, nGroupSize}));
                    condList.add(0);
                    condList.add(new BasicDBObject("$divide", Arrays.asList(new BasicDBObject("$multiply", new Object[] {"$" + MAIN_RESULT_PROJECTION_FIELD + ".he" + g, 100}), new BasicDBObject("$subtract", new Object[] {nGroupSize, "$" + MAIN_RESULT_PROJECTION_FIELD + ".m" + g}))));
                    addFieldsIn.put("hef" + g, new BasicDBObject("$cond", condList)); // heterozygous frequency
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
            projectionFields.put(MAIN_RESULT_PROJECTION_FIELD, new BasicDBObject("$let", let));
        }

        if (fIsMultiRunProject && !annotationMatchList.isEmpty())
            groupFields.put(VariantRunData.SECTION_ADDITIONAL_INFO, new BasicDBObject("$addToSet", "$" + VariantRunData.SECTION_ADDITIONAL_INFO));
            
        if (fIsMultiRunProject)
            pipeline.add(new BasicDBObject("$group", groupFields));
        else if (!fForCounting)
            projectionFields.put("_id", "$_id." + VariantRunDataId.FIELDNAME_VARIANT_ID);
        
        if (!annotationMatchList.isEmpty())
            pipeline.add(new BasicDBObject(new BasicDBObject("$match", new BasicDBObject("$and", annotationMatchList))));  // re-apply this here in case of multiple runs (annotations could be in a different run than the one containing wanted genotypes)

        if (!projectionFields.isEmpty()) {
            projectionFields.put(Assembly.getThreadBoundVariantRefPosPath(), 1);
            projectionFields.put(AbstractVariantData.FIELDNAME_TYPE, 1);
            projectionFields.put(AbstractVariantData.FIELDNAME_KNOWN_ALLELES, 1);
            pipeline.add(new BasicDBObject("$project", projectionFields));
        }
        
        if (addFieldsIn.size() > 0) {
            BasicDBObject addFieldsLet = new BasicDBObject("vars", addFieldsVars);
            addFieldsLet.put("in", addFieldsIn);
            pipeline.add(new BasicDBObject("$addFields", new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD, new BasicDBObject("$let", addFieldsLet))));
        }

        if (fExcludeVariantsWithOnlyMissingData) {  // not all runs are selected: some variants may have no data for the selected samples, we don't want to include them
            ArrayList<BasicDBObject> orList = new ArrayList<BasicDBObject>();
            for (int filteredGroup : filteredGroups) {
                if (req.getMaxMissingData(filteredGroup) == 100)
                    orList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".m" + filteredGroup, new BasicDBObject("$lt", individualToSampleListMap.get(filteredGroup).size())));
            }
            if (!orList.isEmpty())
                finalMatchList.add(new BasicDBObject("$or", orList));
        }
        
        if (finalMatchList.size() > 0)
             pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", finalMatchList)));

        if (nNextCallCount == 1) {
            try { System.err.println(new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString(pipeline)); }
            catch (Exception ignored) {}
        }
        return pipeline;
    }

//    private List<Comparable> buildFullListFromRange(MongoTemplate mongoTemplate, Comparable leftBound, Comparable rightBound) {
//        if (leftBound == null)
//            leftBound = mongoTemplate.findOne(new Query().with(new Sort(Sort.Direction.ASC, "_id")), VariantData.class).getId();
//        if (rightBound == null)
//            rightBound = mongoTemplate.findOne(new Query().with(new Sort(Sort.Direction.DESC, "_id")), VariantData.class).getId();
//
//        ArrayList<Comparable> result = new ArrayList<>();
//        String leftAsString = leftBound.toString(), rightAsString = rightBound.toString();
//        if (ObjectId.isValid(leftAsString) && ObjectId.isValid(rightAsString) && leftAsString.substring(0, 18).equals(rightAsString.substring(0, 18)))
//        {
//            int nCurrentId = Integer.parseInt(leftAsString.substring(18, 24), 16);
//            while (nCurrentId <= Integer.parseInt(rightAsString.substring(18, 24), 16))
//                result.add(new ObjectId(leftAsString.substring(0, 18) + Integer.toHexString(nCurrentId++)));
//        }
//        else
//        {
//            Query q = new Query().with(new Sort(Sort.Direction.ASC, "_id"));
//            q.fields().include("_id");
//            // not finished implementing method (functionality non needed)
//        }
//        return result;
//    }

    public static boolean areObjectIDsConsecutive(ObjectId first, ObjectId second)
    {
//        if (first == null || second == null)
//            return false;

        String firstAsString = first.toHexString(), secondAsString = second.toHexString();
        if (!firstAsString.substring(0, 18).equals(secondAsString.substring(0, 18)))
            return false;
        
        return 1 + Integer.parseInt(firstAsString.substring(18, 24), 16) == Integer.parseInt(secondAsString.substring(18, 24), 16);
    }

//    public static Document tryAndShrinkIdList(String pathToVariantId, Collection<Comparable> idCollection, int nShrinkThreshold)
//    {
//        if (idCollection.size() >= 300000)
//            try
//            {
//        //        long b4 = System.currentTimeMillis();
//        //        SortedSet<Comparable> idSet = SortedSet.class.isAssignableFrom(idCollection.getClass()) ? (SortedSet<Comparable>) idCollection : new TreeSet<Comparable>(idCollection);
//        //        System.out.println("sorting took " + (System.currentTimeMillis() - b4));
//                BasicDBList orList = new BasicDBList();
//                ArrayList<ObjectId> inIdList = new ArrayList<>(), rangeIdList = new ArrayList<>();
//                
//                ObjectId previousId = null;
//                for (Comparable id : idCollection)
//                {
//                    ObjectId currentId = (ObjectId) id;
//                    if (previousId == null || areObjectIDsConsecutive(previousId, currentId))
//                        rangeIdList.add(currentId);
//                    else
//                    {
//                        if (rangeIdList.size() >= nShrinkThreshold)
//                        {    // replace list with a range
//                            BasicDBList chunkMatchAndList = new BasicDBList();
//                            chunkMatchAndList.add(new Document(pathToVariantId, new Document("$gte", rangeIdList.get(0))));
//                            chunkMatchAndList.add(new Document(pathToVariantId, new Document("$lte", rangeIdList.get(rangeIdList.size() - 1))));
//                            orList.add(new Document("$and", chunkMatchAndList));
//                        }
//                        else
//                            inIdList.addAll(rangeIdList);    // range is too small, keep the list
//        
//                        rangeIdList.clear();
//                        rangeIdList.add(currentId);
//                    }
//                    previousId = currentId;
//                }
//                inIdList.addAll(rangeIdList);
//        
//                if (inIdList.size() > 0 || orList.size() == 0)
//                    orList.add(new Document(pathToVariantId, new Document("$in", inIdList)));
//        
//                return orList.size() > 1 ? new Document("$or", orList) : (Document) orList.iterator().next();
//            }
//            catch (ClassCastException cce)
//            {
//                if (!cce.getMessage().contains("ObjectId"))
//                    throw cce;    // otherwise it simply means IDs are of a different type, in which case we can't shrink the collection
//            }
////        else
////        {
////            LOG.debug("Didn't shrink id collection (" + idCollection.size() + " records only)");
////        }
//        
//        return new Document(pathToVariantId, new Document("$in", idCollection));    // not shrinked
//    }
    
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
        for (int i = 0; i < gsvr.getNumberGroups(); i++) {
            if (gsvr.isDiscriminate(i) || !gsvr.getGtPattern(i).equals(GigwaSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL) || gsvr.getMinHeZ(i) > 0 || gsvr.getMaxHeZ(i) < 100 || gsvr.getMinMissingData(i) > 0 || gsvr.getMaxMissingData(i) < 100 || gsvr.getMinMaf(i) > 0 || gsvr.getMaxMaf(i) < 50)
                result.add(i);
        }
        return result;
    }
    
    static public List<Integer> getGroupsForWhichToFilterOnGenotypingOrAnnotationData(GigwaSearchVariantsRequest gsvr, boolean fConsiderFielThresholds)
    {
        List<Integer> result = getGroupsForWhichToFilterOnGenotypingData(gsvr, fConsiderFielThresholds);
 
        if (result.size() == 0 && (gsvr.getGeneName().length() > 0 || gsvr.getVariantEffect().length() > 0))
            result.add(0);    // needed at least for filtering on annotation data or distinguish records according to project id

        /*FIXME: this should also force filtering on VRD in cases where only some runs of the selected project are involved*/
        return result;
    }
}