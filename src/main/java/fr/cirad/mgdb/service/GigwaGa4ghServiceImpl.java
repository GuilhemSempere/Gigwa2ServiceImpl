/*******************************************************************************
 * GIGWA - Service implementation
 * Copyright (C) 2016 - 2025, <CIRAD> <IRD>
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
package fr.cirad.mgdb.service;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.lang.reflect.InvocationTargetException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.ejb.ObjectNotFoundException;
import javax.servlet.ServletContext;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.HttpSession;

import org.apache.avro.AvroRemoteException;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.ga4gh.methods.GAException;
import org.ga4gh.methods.ListReferenceBasesRequest;
import org.ga4gh.methods.ListReferenceBasesResponse;
import org.ga4gh.methods.ReferenceMethods;
import org.ga4gh.methods.SearchCallSetsRequest;
import org.ga4gh.methods.SearchCallSetsResponse;
import org.ga4gh.methods.SearchReferenceSetsRequest;
import org.ga4gh.methods.SearchReferenceSetsResponse;
import org.ga4gh.methods.SearchReferencesRequest;
import org.ga4gh.methods.SearchReferencesResponse;
import org.ga4gh.methods.SearchVariantSetsRequest;
import org.ga4gh.methods.SearchVariantSetsResponse;
import org.ga4gh.methods.SearchVariantsRequest;
import org.ga4gh.methods.VariantMethods;
import org.ga4gh.models.AlleleLocation;
import org.ga4gh.models.AnalysisResult;
import org.ga4gh.models.Call;
import org.ga4gh.models.Call.Builder;
import org.ga4gh.models.CallSet;
import org.ga4gh.models.HGVSAnnotation;
import org.ga4gh.models.OntologyTerm;
import org.ga4gh.models.Reference;
import org.ga4gh.models.ReferenceSet;
import org.ga4gh.models.TranscriptEffect;
import org.ga4gh.models.Variant;
import org.ga4gh.models.VariantAnnotation;
import org.ga4gh.models.VariantSet;
import org.ga4gh.models.VariantSetMetadata;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;
import org.springframework.security.core.Authentication;
import org.springframework.security.core.authority.SimpleGrantedAuthority;
import org.springframework.stereotype.Component;
import org.springframework.util.FileSystemUtils;
import org.springframework.web.context.ServletContextAware;

import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import com.mongodb.MongoCommandException;
import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;
import com.mongodb.client.model.Aggregates;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.individualoriented.AbstractIndividualOrientedExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.importing.SequenceImport;
import fr.cirad.mgdb.importing.VcfImport;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.CachedCount;
import fr.cirad.mgdb.model.mongo.maintypes.DBVCFHeader;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.Sequence;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.Run;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongo.subtypes.VariantRunDataId;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.model.GigwaSearchCallSetsRequest;
import fr.cirad.model.GigwaSearchReferencesRequest;
import fr.cirad.model.GigwaSearchVariantsExportRequest;
import fr.cirad.model.GigwaSearchVariantsResponse;
import fr.cirad.model.MgdbSearchVariantsRequest;
import fr.cirad.security.base.IRoleDefinition;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.AppConfig;
import fr.cirad.tools.ExperimentalFeature;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.SessionAttributeAwareThread;
import fr.cirad.tools.mgdb.GenotypingDataQueryBuilder;
import fr.cirad.tools.mgdb.VariantQueryBuilder;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import fr.cirad.tools.query.GroupedExecutor;
import fr.cirad.tools.query.GroupedExecutor.TaskWrapper;
import fr.cirad.tools.security.base.AbstractTokenManager;
import fr.cirad.utils.Constants;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;

/**
 *
 * Implementation of Gigwa's data-related functionalities
 *
 * GA4GH / Gigwa concept equivalence:
 * GAReferenceSet = Module / Database, GAVariantSet = Project,
 * GAReference = Sequence, GACallSet = Individual, GAVariant = Variant
 *
 * @author adrien petel, guilhem sempere
 */
@Component
public class GigwaGa4ghServiceImpl implements IGigwaService, VariantMethods, ReferenceMethods, ServletContextAware {

    /**
     * logger
     */
    protected static final Logger LOG = Logger.getLogger(GigwaGa4ghServiceImpl.class);

//    /**
//     * The Constant SEQLIST_FOLDER.
//     */
//    static final public String SEQLIST_FOLDER = "selectedSeqs";

    /**
     * The Constant EXPORT_EXPIRATION_DELAY_MILLIS.
     */
    static final private long EXPORT_EXPIRATION_DELAY_MILLIS = 1000 * 60 * 60 * 24 * 2; /* 2 days */
    static final private long DDL_EXPORT_EXPIRATION_DELAY_MILLIS = 1000 * 60 * 60 ; /* 1 hour */

    /**
     * The Constant TMP_OUTPUT_FOLDER.
     */
    static public final String TMP_OUTPUT_FOLDER = "tmpOutput";
    static public final String TMP_OUTPUT_DDL_FOLDER = "ddl_tmpOutput";
	static public final String TMP_OUTPUT_EXTRACTION_FOLDER = "extraction";
	
    /**
     * The Constant FRONTEND_URL.
     */
    static final public String FRONTEND_URL = "genofilt";

    static final protected HashMap<String, String> annotationField = new HashMap<>();


    private Boolean fAllowDiskUse = null;

    @Autowired AbstractTokenManager tokenManager;

    @Autowired private AppConfig appConfig;

    private HashSet<String> hostsNotSupportingMergeOperator = new HashSet<>();

    @Autowired private MgdbDao mgdbDao;

	private ServletContext servletContext;

    /**
     * number format instance
     */
    static protected NumberFormat nf = NumberFormat.getInstance();
    
//    static protected ExecutorService executor;
    
    static {
        nf.setMaximumFractionDigits(4);

        annotationField.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + 1, "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + 1);
        annotationField.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + 2, "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + 2);
        annotationField.put(VariantRunData.SECTION_ADDITIONAL_INFO + "." + 1, "$" + VariantRunData.SECTION_ADDITIONAL_INFO);
    }

    @Override
    public List<String> listVariantTypesSorted(String sModule, int projId) {
        List<String> variantTypesArray = new ArrayList<>(MgdbDao.getVariantTypes(MongoTemplateManager.get(sModule), projId));
        Collections.sort(variantTypesArray, new AlphaNumericComparator<String>());
        return variantTypesArray;
    }

    @Override
    public Collection<String> listModules() {
        return MongoTemplateManager.getAvailableModules();
    }

    @Override
    public int getProjectPloidyLevel(String sModule, int projId) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_PLOIDY_LEVEL);
           q.addCriteria(Criteria.where("_id").is(projId));
        GenotypingProject proj = mongoTemplate.findOne(q, GenotypingProject.class);
        return proj.getPloidyLevel();
    }

    @Override
    public TreeSet<String> searchableAnnotationFields(String sModule, int projId) {
        /* This may be more efficient by looking at the VCF header instead */
        TreeSet<String> result = new TreeSet<>();
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Query q = new Query(Criteria.where("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID).is(projId));
        q.limit(3);
        q.fields().include(VariantRunData.FIELDNAME_SAMPLEGENOTYPES);
        Iterator<VariantRunData> it = mongoTemplate.find(q, VariantRunData.class).iterator();
        while (it.hasNext())
        {
            VariantRunData vrd = it.next();
            Collection<SampleGenotype> spGTs = vrd.getSampleGenotypes().size() > 100 ? new ArrayList<>(vrd.getSampleGenotypes().values()).subList(0,  100) : vrd.getSampleGenotypes().values();
            for (HashMap<String, Object> annotationMap : spGTs.stream().map(sg -> sg.getAdditionalInfo()).collect(Collectors.toList()))
                for (String aiKey : annotationMap.keySet())
                    if (Number.class.isAssignableFrom(annotationMap.get(aiKey).getClass()))
                         result.add(aiKey);
        }
        return result;
    }

    @Override
    public TreeSet<String> getProjectEffectAnnotations(String sModule, int projId) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_EFFECT_ANNOTATIONS);
           q.addCriteria(Criteria.where("_id").is(projId));
        GenotypingProject proj = mongoTemplate.findOne(q, GenotypingProject.class);
        return proj.getEffectAnnotations();
    }

    @Override
    public Collection<Integer> getDistinctAlleleCounts(String sModule, Integer projId) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        return mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(GenotypingProject.class)).distinct(GenotypingProject.FIELDNAME_ALLELE_COUNTS, projId == null ? null : new BasicDBObject("_id", projId), Integer.class).into(new ArrayList<>());
    }

    @Override
    public Map<Integer, String> getProjectIdToNameMap(String sModule) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Map<Integer, String> result = new LinkedHashMap<>();
        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_NAME);
        for (GenotypingProject proj : mongoTemplate.find(q, GenotypingProject.class)) {
            result.put(proj.getId(), proj.getName());
        }
        return result;
    }

    @Override
    public String getQueryForGenotypePattern(String gtPattern) {
        return GenotypingDataQueryBuilder.getGenotypePatternToQueryMap().get(gtPattern);
    }

    @Override
    public List<String> listIndividualsInAlphaNumericOrder(String sModule, int project) {
        List<String> indArray = null;
        try {
            indArray = new ArrayList(MgdbDao.getProjectIndividuals(sModule, project));
        } catch (ObjectNotFoundException ex) {
            java.util.logging.Logger.getLogger(GigwaGa4ghServiceImpl.class.getName()).log(Level.SEVERE, null, ex);
        }
        Collections.sort(indArray, new AlphaNumericComparator());
        return indArray;
    }

    @Override
    public List<String> listSequences(HttpServletRequest request, String sModule, int projId) {
        List<String> result = new ArrayList<String>(MongoTemplateManager.get(sModule).findById(projId, GenotypingProject.class).getContigs(Assembly.getThreadBoundAssembly()));

//        List<String> externallySelectedSequences = getSequenceIDsBeingFilteredOn(request.getSession(), sModule);
//        /* first try to use a list that may have been defined on in a different section of the application (although it may not be limited to the given project) */
//        if (externallySelectedSequences != null) {
//            result = (List<String>) CollectionUtils.intersection(result, externallySelectedSequences);
//        }

        if (result != null) {
            Collections.sort(result, new AlphaNumericComparator());
        }
        return result;
    }

    private String isSearchedDatasetReasonablySized(MgdbSearchVariantsRequest gsvr) throws Exception
    {
        List<Integer> groupsForWhichToFilterOnGenotypingData = VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingData(gsvr, false);
        if (groupsForWhichToFilterOnGenotypingData.isEmpty())
            return null;    // no genotyping data filtering involved

        String info[] = Helper.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

        List<List<String>> callsetIds = gsvr.getAllCallSetIds();

        int nIndCount = 0;
        for (int i = 0; i < callsetIds.size(); i++)
        	if (groupsForWhichToFilterOnGenotypingData.contains(i))
        		nIndCount += callsetIds.get(i).size();
        if (nIndCount == 0)
        	nIndCount = MgdbDao.getProjectIndividuals(info[0], projId).size();

        int nMaxBillionGenotypesInvolved = 1;    // default
        try
        {
            nMaxBillionGenotypesInvolved = Integer.parseInt(appConfig.get("maxSearchableBillionGenotypes"));
        }
        catch (Exception ignored)
        {}

        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Query q = new Query();
        q.fields().include(Assembly.getThreadBoundProjectContigsPath());
        q.addCriteria(Criteria.where("_id").is(projId));
        GenotypingProject proj = mongoTemplate.findOne(q, GenotypingProject.class);

        int nSelectedSeqCount = gsvr.getReferenceName() == null || gsvr.getReferenceName().length() == 0 ? proj.getContigs(Assembly.getThreadBoundAssembly()).size() : gsvr.getReferenceName().split(";").length;
        if (nSelectedSeqCount == 1)
            return null;    // we can't expect user to select less than a single sequence

        int nAvgVariantsPerSeq = (int) (Helper.estimDocCount(mongoTemplate, VariantData.class) / Math.max(1, proj.getContigs(Assembly.getThreadBoundAssembly()).size()));
        BigInteger maxSeqCount = BigInteger.valueOf(1000000000).multiply(BigInteger.valueOf(nMaxBillionGenotypesInvolved)).divide(BigInteger.valueOf(nAvgVariantsPerSeq).multiply(BigInteger.valueOf(nIndCount)));
        int nMaxSeqCount = Math.max(1, maxSeqCount.intValue());

        if (nSelectedSeqCount <= nMaxSeqCount)
            return null;

        return "This is a big database. Given the number of selected individuals, you may only work on a maximum of " + nMaxSeqCount + " sequence(s) at a time.";
    }

    private void applyPreFiltering(List<BasicDBObject> genotypingDataPipeline, boolean fMongoOnSameServer, MongoCollection<Document> varColl) {
        BasicDBObject initialMatch = (BasicDBObject) genotypingDataPipeline.get(0).get("$match");
        if (initialMatch != null) {    // initialMatchForVariantColl will be the one applied to variants collection when pre-filtering
            BasicDBList initialMatchForVariantColl = (BasicDBList) ((BasicDBList) initialMatch.get("$and")).clone();
            boolean fMultiProjectDB = false;
            if (initialMatchForVariantColl != null) {
                List<DBObject> toAdd = new ArrayList<>(), toRemove = new ArrayList<>();
                for (Object filter : initialMatchForVariantColl) {
                    Object variantIdFilter = ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID);
                    if (variantIdFilter != null) {
                        toAdd.add(new BasicDBObject("_id", variantIdFilter));
                        toRemove.add((DBObject) filter);
                    }
                    else if (null != ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID)) {
                        toRemove.add((DBObject) filter);    // no project info to filter on in the variants collection
                        fMultiProjectDB = true;
                    }
                    else {
                    	Object runFilter = ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_RUNNAME);
                    	if (runFilter != null) {
	                    	toAdd.add(new BasicDBObject(VariantData.FIELDNAME_RUNS + "." + Run.FIELDNAME_RUNNAME, runFilter));
	                        toRemove.add((DBObject) filter);    // no project info to filter on in the variants collection
                        }
                    }
                }
                initialMatchForVariantColl.addAll(toAdd);
                initialMatchForVariantColl.removeAll(toRemove);
            }
            
            if (fMongoOnSameServer) {    // always worth pre-filtering
                MongoCursor<Document> variantCursor = varColl.find(new BasicDBObject("$and", initialMatchForVariantColl)).projection(new BasicDBObject("_id", 1)).iterator();
                List<Comparable> chunkPreFilteredIDs = new ArrayList<>();
                while (variantCursor.hasNext())
                	chunkPreFilteredIDs.add((Comparable) variantCursor.next().get("_id"));
                if (chunkPreFilteredIDs.size() == 0)
                	genotypingDataPipeline.clear();    // no variants match indexed part of the query: skip chunk
                else {    // DB server is the same machine as web server: $in operator will not be expensive
                    if (!fMultiProjectDB)    // for single project dbs, $in is equivalent to original query, otherwise only a pre-filter
                    	genotypingDataPipeline.remove(0);
                    genotypingDataPipeline.add(0, new BasicDBObject("$match", new BasicDBObject("$and", new BasicDBList() {{ add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", chunkPreFilteredIDs)));}} )));
                }
            }
            else if (varColl.countDocuments(new BasicDBObject("$and", initialMatchForVariantColl)) == 0)    // only try and use pre-filtering to avoid executing genotyping data queries on irrelevant chunks
               	genotypingDataPipeline.clear();    // no variants match indexed part of the query: skip chunk
        }
    }
    
    @Override
    public long countVariants(MgdbSearchVariantsRequest gsvr, boolean fSelectionAlreadyExists) throws Exception {
        String info[] = Helper.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

        boolean fGotTokenManager = tokenManager != null;    // if null, we are probably being invoked via unit-test
        String token = !fGotTokenManager ? Helper.convertToMD5(String.valueOf(System.currentTimeMillis())) /* create a mock up token */ : tokenManager.readToken(gsvr.getRequest());

        ProgressIndicator progress = ProgressIndicator.get(token);    // it may already exist (if we're being called by findVariants for example)
        if (progress == null) {
            progress = new ProgressIndicator(token, new String[0]);
            ProgressIndicator.registerProgressIndicator(progress);
        }
        String sizeProblemMsg = gsvr.shallApplyMatrixSizeLimit() ? isSearchedDatasetReasonablySized(gsvr) : null;
        if (sizeProblemMsg != null)
        {
            progress.setError(sizeProblemMsg);
            return 0;
        }

        String queryKey = getQueryKey(gsvr);
        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        
        Long count = CachedCount.getCachedCount(mongoTemplate, queryKey, null);
        LOG.debug((count == null ? "new" : "existing") + " queryKey hash: " + queryKey);
        if (count == null)
        {
            long before = System.currentTimeMillis();
            progress.addStep("Counting matching variants");

            List<String> alleleCountList = gsvr.getAlleleCount().length() == 0 ? null : Arrays.asList(gsvr.getAlleleCount().split(";"));

            GenotypingProject genotypingProject = mongoTemplate.findById(projId, GenotypingProject.class);
            if (genotypingProject.getAlleleCounts().size() != 1 || genotypingProject.getAlleleCounts().iterator().next() != 2) {    // Project does not only have bi-allelic data: make sure we can apply MAF filter on selection
                boolean fExactlyOneNumberOfAllelesSelected = alleleCountList != null && alleleCountList.size() == 1;
                boolean fBiAllelicSelected = fExactlyOneNumberOfAllelesSelected && "2".equals(alleleCountList.get(0));
                for (int i = 0; i < gsvr.getNumberGroups(); i++) 
	                if (!fBiAllelicSelected && (gsvr.getMaxMaf(i) < 50 || gsvr.getMinMaf(i) > 0)) {
	                    progress.setError("MAF is only supported on biallelic data!");
	                    return 0l;
                }
            }

            MongoCollection<Document> varColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class));
            List<Integer> filteredGroups = VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gsvr, false);

            VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gsvr/*, !fGotTokenManager ? null : getSequenceIDsBeingFilteredOn(gsvr.getRequest().getSession(), sModule)*/, false);
        	Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
            if (filteredGroups.isEmpty()) { // filtering on variant features only: we just need a count
                if (variantDataQueries.isEmpty())
                    count = mongoTemplate.count(new Query(), VariantData.class);    // no filter whatsoever
                else {
                    count = 0l;
                    count = variantDataQueries.parallelStream()
                            .map(req -> varColl.countDocuments(req.isEmpty() ? new BasicDBObject() : new BasicDBObject("$and", req)))
                            .reduce(count, (accumulator, _item) -> accumulator + _item);
                }
            }

            if (count != null)
            	mongoTemplate.save(new CachedCount(queryKey, Arrays.asList(count)));
            else
            {    // filter on genotyping data
                boolean fPreFilterOnVarColl = false, fMongoOnSameServer = MongoTemplateManager.isModuleOnLocalHost(sModule);

                //in this case, there is only one variantQueryDBList (no filtering on variant ids)
                BasicDBList variantQueryDBList = !variantDataQueries.isEmpty() ? variantDataQueries.iterator().next() : new BasicDBList();

                if (variantQueryDBList.size() > 0)
                {
                    Number avgObjSize = (Number) mongoTemplate.getDb().runCommand(new BasicDBObject("collStats", mongoTemplate.getCollectionName(VariantRunData.class))).get("avgObjSize");
                    if (avgObjSize.doubleValue() >= 10240)
                    {    // it may be worth pre-filtering data on variant collection because filtering speed on the run collection is affected by the document size
                        long totalCount = mongoTemplate.count(new Query(), VariantData.class), preFilterCount = varColl.countDocuments(new BasicDBObject("$and", variantQueryDBList));
                        fPreFilterOnVarColl = preFilterCount <= totalCount*(fMongoOnSameServer ? .85 : .45);    // only pre-filter if less than a given portion of the total variants are to be retained
                        if (fPreFilterOnVarColl)
                        	LOG.debug("Pre-filtering data on variant collection");
                    }
                }

                GenotypingDataQueryBuilder genotypingDataQueryBuilder = new GenotypingDataQueryBuilder(gsvr, variantQueryDBList, true);
                final int nChunkCount = genotypingDataQueryBuilder.getNumberOfQueries();
                final List<Integer> shuffledChunkIndexes = genotypingDataQueryBuilder.shuffleChunkOrder();
                try
                {
                    if (nChunkCount > 1)
                        LOG.debug("Query split into " + nChunkCount);

                    final Long[] partialCountArray = new Long[nChunkCount];
                    final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();
                    final AtomicInteger finishedThreadCount = new AtomicInteger(0);

                    int i = -1;
                    String taskGroup = "count_" + System.currentTimeMillis() + "_" + token;
//                    Long b4 = System.currentTimeMillis();

                    ExecutorService executor = MongoTemplateManager.getExecutor(sModule);
                    Integer nextConcurrentThreadCountReevaluationChunk = executor instanceof GroupedExecutor ? null : ((ThreadPoolExecutor) executor).getCorePoolSize();

                    while (genotypingDataQueryBuilder.hasNext()) {
                    	List<BasicDBObject> genotypingDataPipeline = genotypingDataQueryBuilder.next();
                        final int chunkIndex = shuffledChunkIndexes.get(++i);

                        // Now the $group operation, used for counting
                        genotypingDataPipeline.add(new BasicDBObject("$count", "count"));

                        if (progress.isAborted())
                            return 0l;

                        final ProgressIndicator finalProgress = progress;
                        Thread queryThread = new Thread() {
                            @Override
                            public void run() {
                            	if (finalProgress.getError() == null && !finalProgress.isAborted())
	                                try {
                                		applyPreFiltering(genotypingDataPipeline, fMongoOnSameServer, varColl);
                                		if (genotypingDataPipeline.isEmpty()) {
                                            partialCountArray[chunkIndex] = 0l;	// no variants match indexed part of the query: skip chunk
                                            finalProgress.setCurrentStepProgress((short) (finishedThreadCount.incrementAndGet() * 100 / nChunkCount));
                                            return;
                                		}
	                                    MongoCursor<Document> it = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(genotypingDataPipeline).allowDiskUse(true).iterator();
	                                    partialCountArray[chunkIndex] = it.hasNext() ? ((Number) it.next().get("count")).longValue() : 0;
	                                    finalProgress.setCurrentStepProgress((short) (finishedThreadCount.incrementAndGet() * 100 / nChunkCount));
	                                    it.close();
	                                }
	                                catch (Throwable t) {
	                                	if (finalProgress.getError() == null && !finalProgress.isAborted()) {
		                                    LOG.error("Error counting variants", t);
		                                    finalProgress.setError(t.getMessage());
	                                	}
	                                }
                            	finally {
                                    genotypingDataPipeline.clear();    // release memory (VERY IMPORTANT)
                            	}
                            	else {
                            		for (Future<Void> f : threadsToWaitFor)
                            			f.cancel(true);
                            		LOG.debug("Cancelled query threads for process " + finalProgress.getProcessId());
                            	}
                            }
                        };
                        
                        Future<Void> queryRunnable = (Future<Void>) executor.submit(new TaskWrapper(taskGroup, queryThread));
                        threadsToWaitFor.add(queryRunnable);
                        
                        if (nextConcurrentThreadCountReevaluationChunk != null) {
                        	ThreadPoolExecutor threadExecutor = ((ThreadPoolExecutor) executor);
                            int nNConcurrentThreads = threadExecutor.getCorePoolSize();
                            if (threadsToWaitFor.size() == nextConcurrentThreadCountReevaluationChunk) {
	                        	queryRunnable.get();

		                        // regulate number of concurrent threads
	                            int nRunningThreadCount = threadExecutor.getActiveCount();
	                            if (nRunningThreadCount > nNConcurrentThreads * .5)
	                                nNConcurrentThreads /= 1.5;
	                            else if (nRunningThreadCount < nNConcurrentThreads * .25)
	                                nNConcurrentThreads *= 1.5;
	                            nNConcurrentThreads = Math.min(MongoTemplateManager.MAXIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS, Math.max(MongoTemplateManager.MINIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS, nNConcurrentThreads));
	                            threadExecutor.setCorePoolSize(nNConcurrentThreads);
	                            nextConcurrentThreadCountReevaluationChunk = threadsToWaitFor.size() + nNConcurrentThreads;
//	                            System.out.println(nRunningThreadCount + " / " + threadsToWaitFor.size() + " -> " + nNConcurrentThreads);
                        	}
                        }
                    }
                    
	                if (nextConcurrentThreadCountReevaluationChunk == null) {
//		                LOG.debug("Submitting queries on " + ((GroupedExecutor) executor).getCorePoolSize() + " threads took " + (System.currentTimeMillis() - b4) + "ms");
	                	((GroupedExecutor) executor).shutdown(taskGroup);	// important to be sure that all tasks in the group are executed before the queue purges it
	                }

                    for (Future<Void> t : threadsToWaitFor) // wait for all threads before moving to next phase
                    	t.get();
                    
                    count = 0l;
                    if (progress.getError() == null && !progress.isAborted()) {
	                    progress.setCurrentStepProgress(100);
	
	                    for (Long partialCount : partialCountArray)
	                    	count += partialCount;
	
	                    mongoTemplate.save(new CachedCount(queryKey, Arrays.asList(partialCountArray)));
                    }
                    if (nextConcurrentThreadCountReevaluationChunk != null)
	                	executor.shutdownNow();
                }
                catch (CancellationException ignored) {}
                catch (Exception e) {
                	if (progress.getError() == null && !progress.isAborted())
                		progress.setError(e.getMessage());
                    LOG.warn("Error searching variants", e);
                }
            }
            LOG.info("countVariants found " + count + " results in " + (System.currentTimeMillis() - before) / 1000d + "s");
        }

        if (progress.getError() != null || progress.isAborted())
            return 0l;

        progress.markAsComplete();
        return count;
    }

//	/**
//     * Gets the temporary variant collection.
//     *
//     * @param sModule the module
//     * @param processID the process id
//     * @param fEmptyItBeforeHand whether or not to empty it beforehand
//     * @return the temporary variant collection
//	 * @throws InterruptedException 
//     */
//    public MongoCollection<Document> getTemporaryVariantCollection(String sModule, String processID, boolean fEmptyItBeforeHand) throws InterruptedException {
//        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
//        MongoCollection<Document> tmpColl = mongoTemplate.getCollection(MongoTemplateManager.TEMP_COLL_PREFIX + Helper.convertToMD5(processID));
//        if (fEmptyItBeforeHand) {
//
////            ArrayList<StackTraceElement> keptStackTraceElements = new ArrayList<>();
////            Exception e = new Exception("Check stack trace");
////            for (StackTraceElement ste : e.getStackTrace())
////                if (ste.toString().startsWith("fr.cirad."))
////                    keptStackTraceElements.add(ste);
////            e.setStackTrace(keptStackTraceElements.toArray(new StackTraceElement[keptStackTraceElements.size()]));
////            LOG.debug("Dropping " + sModule + "." + tmpColl.getName() + " from getTemporaryVariantCollection", e);
//
//            tmpColl.drop();
//            MgdbDao.ensurePositionIndexes(mongoTemplate, Arrays.asList(tmpColl), false, false);    // make sure we have indexes defined as required in v2.4
//        }
//        return tmpColl;
//    }

    @Override
    public long findVariants(MgdbSearchVariantsRequest gsvr) throws Exception {
        String token = tokenManager.readToken(gsvr.getRequest());

        final ProgressIndicator progress = new ProgressIndicator(token, new String[0]);
        ProgressIndicator.registerProgressIndicator(progress);
        String sizeProblemMsg = gsvr.shallApplyMatrixSizeLimit() ? isSearchedDatasetReasonablySized(gsvr) : null;
        if (sizeProblemMsg != null) {
            progress.setError(sizeProblemMsg);
            return 0;
        }

        progress.addStep("Finding matching variants");

        String info[] = Helper.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];
        String queryKey = getQueryKey(gsvr);

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        String sMongoHost = MongoTemplateManager.getModuleHost(sModule);

        MongoCollection<Document> cachedCountCollection = mongoTemplate.getCollection(mongoTemplate.getCollectionName(CachedCount.class));
        MongoCursor<Document> countCursor = cachedCountCollection.find(new BasicDBObject("_id", queryKey)).iterator();

        final Object[] partialCountArray = !countCursor.hasNext() ? null : ((List<Object>) countCursor.next().get(MgdbDao.FIELD_NAME_CACHED_COUNT_VALUE)).toArray();
        final LinkedHashMap<Integer, Long> partialCountMap = new LinkedHashMap<>(); // progress display will be more accurate if we skip empty chunks
        long nTotalCount = 0;
        if (partialCountArray != null) {
            for (int i=0; i<partialCountArray.length; i++) {
                long n = (long) partialCountArray[i];
                if (n != 0) {
                    partialCountMap.put(i, n);
                    nTotalCount += n;
                }
            }
            if (nTotalCount == 0) {
                progress.markAsComplete();
                return 0;
            }
        }
        LOG.debug((partialCountArray == null ? "new" : "existing") + " queryKey hash: " + queryKey);

        List<Integer> filteredGroups = VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gsvr, false);
        boolean fNeedToCreateTempColl = !filteredGroups.isEmpty() 
                || (gsvr.getSelectedVariantIds() != null && !gsvr.getSelectedVariantIds().isEmpty());

        long before = System.currentTimeMillis();

        if (fNeedToCreateTempColl) {
            final MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(sModule, progress.getProcessId(), true, false, false);

            VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gsvr/*, getSequenceIDsBeingFilteredOn(gsvr.getRequest().getSession(), sModule)*/, false);
            Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
            
	        if (filteredGroups.isEmpty()) {	// filtering on Variant collection (filtering on variant IDs or run)
	            MongoCollection<Document> varColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class));
	            final AtomicInteger nProcessedChunkCount = new AtomicInteger(0);
	
	            variantDataQueries.parallelStream().forEach(req -> {
	                varColl.aggregate(Arrays.asList(
	                    new BasicDBObject("$match", new BasicDBObject("$and", req)),
	                    new BasicDBObject("$merge", new BasicDBObject("into", tmpVarColl.getNamespace().getCollectionName()).append("whenMatched", "fail"))
	                )).toCollection();
	
	                progress.setCurrentStepProgress(nProcessedChunkCount.incrementAndGet() * 100 / variantDataQueries.size());
	            });
	        }
	        else {   // filter on genotyping data
	            final ArrayList<Future<Void>> threadsToWaitFor = new ArrayList<>();
	            final AtomicInteger finishedThreadCount = new AtomicInteger(0);
                    
	            Collection<BasicDBList> variantRunDataQueries = varQueryWrapper.getVariantRunDataQueries();
	            
	            //in this case, there is only one variantQueryDBList (no filtering on variant ids)
	            BasicDBList variantQueryDBList = !variantRunDataQueries.isEmpty() ? variantRunDataQueries.iterator().next() : new BasicDBList();
	
	            final GenotypingDataQueryBuilder genotypingDataQueryBuilder = new GenotypingDataQueryBuilder(gsvr, variantQueryDBList, false);
	
	            try {
		                final int nChunkCount = genotypingDataQueryBuilder.getNumberOfQueries();
		                final List<Integer> shuffledChunkIndexes = genotypingDataQueryBuilder.shuffleChunkOrder();
		
		                final Long[] partialCountArrayToFill = partialCountArray == null ? new Long[nChunkCount] : null;
		                if (partialCountArrayToFill != null)
		                    LOG.info("Find without prior count: do both at once");
		                else if (nChunkCount != partialCountArray.length) {
		                    progress.setError("Different number of chunks between counting and listing variant rows!");
		                    LOG.error(progress.getError());
		                    return 0;
		                }
		
		                MongoCollection<Document> varColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class));
		
		                boolean fPreFilterOnVarColl = false, fMongoOnSameServer = MongoTemplateManager.isModuleOnLocalHost(sModule);
		                if (variantQueryDBList.size() > 0) {
		                        Number avgObjSize = (Number) mongoTemplate.getDb().runCommand(new BasicDBObject("collStats", mongoTemplate.getCollectionName(VariantRunData.class))).get("avgObjSize");
		                    if (avgObjSize.doubleValue() >= 10240) {   // it may be worth pre-filtering data on variant collection because filtering speed on the run collection is affected by the document size
		                        long totalCount = mongoTemplate.count(new Query(), VariantData.class), preFilterCount = varColl.countDocuments(new BasicDBObject("$and", variantQueryDBList));
		                        fPreFilterOnVarColl = preFilterCount <= totalCount*(fMongoOnSameServer ? .85 : .45);    // only pre-filter if less than a given portion of the total variants are to be retained
		                        if (fPreFilterOnVarColl)
		                            LOG.debug("Pre-filtering data on variant collection");
		                    }
		                }
		
		                if (nChunkCount > 1)
		                    LOG.debug("Query split into " + nChunkCount);

		                int i = -1;
		                final MongoCollection<Document> vrdColl = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class));
		                final HashMap<Integer, BasicDBList> rangesToCount = partialCountArrayToFill != null ? new HashMap<>() : null;
                        String taskGroup = "find_" + System.currentTimeMillis() + "_" + token;
//                        Long b4 = System.currentTimeMillis();
                        
                        ExecutorService executor = MongoTemplateManager.getExecutor(sModule);
                        Integer nextConcurrentThreadCountReevaluationChunk = executor instanceof GroupedExecutor ? null : ((ThreadPoolExecutor) executor).getCorePoolSize();
                        
		                while (genotypingDataQueryBuilder.hasNext()) {
		                    List<BasicDBObject> genotypingDataPipeline = genotypingDataQueryBuilder.next();
		                    if (progress.isAborted() || progress.getError() != null)
		                        return 0;
		
		                    final int chunkIndex = shuffledChunkIndexes.get(++i);
		                    if (partialCountMap.size() > 0 && !partialCountMap.containsKey(chunkIndex))
		                        continue;	// we know there are no matches in this chunk
		                    		
	                        if (rangesToCount != null) {	// we need to keep track of the count after searching so prepare queries to run on tmp coll
			                    BasicDBObject initialMatch = (BasicDBObject) genotypingDataPipeline.get(0).get("$match");
			                    BasicDBList initialMatchForVariantColl = (BasicDBList) ((BasicDBList) initialMatch.get("$and")).clone();
			                    List<DBObject> toAdd = new ArrayList<>(), toRemove = new ArrayList<>();
			                    for (Object filter : initialMatchForVariantColl) {
			                        Object variantIdFilter = ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID);
			                        if (variantIdFilter != null) {
			                            toAdd.add(new BasicDBObject("_id", variantIdFilter));
			                            toRemove.add((DBObject) filter);
			                        }
			                        else if (null != ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID) || null != ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_RUNNAME)) {
			                            toRemove.add((DBObject) filter);    // no project / run info to filter on in the variants collection
			                        }
			                    }
			                    initialMatchForVariantColl.addAll(toAdd);
			                    initialMatchForVariantColl.removeAll(toRemove);
			                    rangesToCount.put(chunkIndex, initialMatchForVariantColl);
	                        }
		
		                    if (partialCountArray != null)
		                        genotypingDataPipeline.add(new BasicDBObject("$limit", partialCountArray[chunkIndex]));
		                    genotypingDataPipeline.add(new BasicDBObject("$project", new BasicDBObject(VariantData.FIELDNAME_KNOWN_ALLELES, 1).append(Assembly.getThreadBoundVariantRefPosPath(), 1).append(VariantData.FIELDNAME_TYPE, 1)));

	                        Thread queryThread = new Thread() {
	                            @Override
	                            public void run() {
	                            	if (progress.getError() == null && !progress.isAborted()) {
                                		applyPreFiltering(genotypingDataPipeline, fMongoOnSameServer, varColl);
                                		if (genotypingDataPipeline.isEmpty()) {
                                			partialCountArrayToFill[chunkIndex] = 0l;	// no variants match indexed part of the query: skip chunk
                                			progress.setCurrentStepProgress((short) (finishedThreadCount.incrementAndGet() * 100 / (partialCountMap.isEmpty() ? nChunkCount : partialCountMap.size())));
                                            return;
                                		}

	                                    boolean fMergeFailedOnThisChunk = false;
	                                    if (!hostsNotSupportingMergeOperator.contains(sMongoHost))
	                                        try {
	                                            genotypingDataPipeline.add(new BasicDBObject("$merge", new BasicDBObject("into", tmpVarColl.getNamespace().getCollectionName()).append("whenMatched", "fail" /* important (fastest option)*/)));
	                                            vrdColl.aggregate(genotypingDataPipeline).allowDiskUse(true).toCollection();
	                                        }
	                                        catch (Throwable t) {
	                                            if (t instanceof MongoCommandException && t.getMessage().contains("$merge")) {
	                                                hostsNotSupportingMergeOperator.add(sMongoHost);
	                                                fMergeFailedOnThisChunk = true;
	                                                LOG.warn("Disabling use of $merge in creating temporary collections on host " + sMongoHost + " (operator not supported by MongoDB server version)");
	                                            }
	                                            else {
	                                            	if (progress.getError() == null && !progress.isAborted()) {
		                                                LOG.error("Error searching variants", t);
		                                                progress.setError(t.getMessage());
	                                            	}
	                                                return;
	                                            }
	                                        }
		                                if (hostsNotSupportingMergeOperator.contains(sMongoHost))
		                                    try {
	                                            if (fMergeFailedOnThisChunk)
	                                            	genotypingDataPipeline.remove(genotypingDataPipeline.size() - 1);    // remove the $merge step we added
	                                            MongoCursor<Document> genotypingDataCursor = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(genotypingDataPipeline).allowDiskUse(true).iterator();
	                                            final ArrayList<Document> variantsThatPassedRunFilterForThisChunk = new ArrayList<>();
	                                            while (genotypingDataCursor.hasNext())
	                                                variantsThatPassedRunFilterForThisChunk.add(genotypingDataCursor.next());
	
	                                            if (partialCountArrayToFill != null)
	                                                partialCountArrayToFill[chunkIndex] = (long) variantsThatPassedRunFilterForThisChunk.size();
	                                            if (variantsThatPassedRunFilterForThisChunk.size() > 0)
	                                                tmpVarColl.insertMany(variantsThatPassedRunFilterForThisChunk);
	
	                                            genotypingDataCursor.close();
		                                    }
		                                    catch (Exception e) {
		                                    	if (progress.getError() == null && !progress.isAborted()) {
			                                        LOG.error("Error searching variants", e);
			                                        progress.setError(e.getMessage());
		                                    	}
		                                    }
	                                    genotypingDataPipeline.clear();    // release memory
	                                    progress.setCurrentStepProgress((short) (finishedThreadCount.incrementAndGet() * 100 / (partialCountMap.isEmpty() ? nChunkCount : partialCountMap.size())));
//	                                    System.err.println(finishedThreadCount.get() + " / " + progress.getCurrentStepProgress() + " / " + chunkIndex);
	                            	}
	                            	else {
	                            		for (Future<Void> f : threadsToWaitFor)
	                            			f.cancel(true);
	                            		LOG.debug("Cancelled query threads for process " + progress.getProcessId());
	                            	}
	                            }
	                        };

	                        Future<Void> queryRunnable = (Future<Void>) executor.submit(new TaskWrapper(taskGroup, queryThread));
	                        threadsToWaitFor.add(queryRunnable);
	                        
	                        if (nextConcurrentThreadCountReevaluationChunk != null) {
	                        	ThreadPoolExecutor threadExecutor = ((ThreadPoolExecutor) executor);
	                            int nNConcurrentThreads = threadExecutor.getCorePoolSize();
	                            if (threadsToWaitFor.size() == nextConcurrentThreadCountReevaluationChunk) {
		                        	queryRunnable.get();

			                        // regulate number of concurrent threads
		                            int nRunningThreadCount = threadExecutor.getActiveCount();
		                            if (nRunningThreadCount > nNConcurrentThreads * .5)
		                                nNConcurrentThreads /= 1.5;
		                            else if (nRunningThreadCount < nNConcurrentThreads * .25)
		                                nNConcurrentThreads *= 1.5;
		                            nNConcurrentThreads = Math.min(MongoTemplateManager.MAXIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS, Math.max(MongoTemplateManager.MINIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS, nNConcurrentThreads));
		                            threadExecutor.setCorePoolSize(nNConcurrentThreads);
		                            nextConcurrentThreadCountReevaluationChunk = threadsToWaitFor.size() + nNConcurrentThreads;
//		                            System.out.println(nRunningThreadCount + " / " + threadsToWaitFor.size() + " -> " + nNConcurrentThreads);
	                        	}
	                        }
		                }

		                if (nextConcurrentThreadCountReevaluationChunk == null) {
//			                LOG.debug("Submitting queries on " + ((GroupedExecutor) executor).getCorePoolSize() + " threads took " + (System.currentTimeMillis() - b4) + "ms");
		                	((GroupedExecutor) executor).shutdown(taskGroup);	// important to be sure that all tasks in the group are executed before the queue purges it
		                }
		                
	                    for (Future<Void> t : threadsToWaitFor) // wait for all threads before moving to next phase
	                    	t.get();
	
	                    if (progress.getError() == null && !progress.isAborted()) {
	                        progress.setCurrentStepProgress(100);
	                        
	                        if (nextConcurrentThreadCountReevaluationChunk != null)
	                        	((ThreadPoolExecutor) executor).setCorePoolSize(MongoTemplateManager.MINIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS);
	
	                        if (partialCountArrayToFill != null) {    // we don't have a count cache for this query: let's create it
	                            if (!hostsNotSupportingMergeOperator.contains(sMongoHost)) {
	                            	threadsToWaitFor.clear();
	                                taskGroup = "find_count_" + System.currentTimeMillis() + "_" + token;
	                                for (Integer j : rangesToCount.keySet()) {
	        	                        Thread countThread = new Thread() {
	        	                            @Override
	        	                            public void run() {
	                                            partialCountArrayToFill[j] = tmpVarColl.countDocuments(new BasicDBObject("$and", rangesToCount.get(j)));
	                                        }
	                                    };
	                                    threadsToWaitFor.add((Future<Void>) executor.submit(new TaskWrapper(taskGroup, countThread)));
	                                }
	                                
	                                if (nextConcurrentThreadCountReevaluationChunk == null)
	                                	((GroupedExecutor) executor).shutdown(taskGroup);	// important to be sure that all tasks in the group are executed before the queue purges it

	                                for (Future<Void> t : threadsToWaitFor) // wait for all threads before moving to next phase
	                                	t.get();
	                            }
	                            mongoTemplate.save(new CachedCount(queryKey, Arrays.asList(partialCountArrayToFill)));
	                        }
	                    }
                        if (nextConcurrentThreadCountReevaluationChunk != null)
    	                	executor.shutdownNow();
	                }
                    catch (CancellationException ignored) {}
                    catch (Exception e) {
                        LOG.debug("Error searching variants", e);
                    }
	        	}
        }
        else
        	dropTempColl(sModule, progress.getProcessId());	// important because if it remains then a subsequent export will use that instead of working directly on runs

        if (progress.isAborted() || progress.getError() != null)
            return 0;
        
        if (partialCountArray == null)
            nTotalCount = countVariants(gsvr, true);
        LOG.info("findVariants found " + nTotalCount + " results in " + (System.currentTimeMillis() - before) / 1000d + "s");

        progress.markAsComplete();
        return nTotalCount;
    }

    @Override
    public void exportVariants(GigwaSearchVariantsExportRequest gsver, String token, HttpServletResponse response) throws Exception
    {
        String processId = "export_" + token;
        final ProgressIndicator progress = new ProgressIndicator(processId, new String[0]);
        ProgressIndicator.registerProgressIndicator(progress);

//        new Thread() { public void run() {
//                try {
//                    cleanupExpiredExportData(gsver.getRequest().servletContext);
//                } catch (IOException e) {
//                    LOG.error("Unable to cleanup expired export files", e);
//                }
//            }
//        }.start();

        String info[] = Helper.getInfoFromId(gsver.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

        long before = System.currentTimeMillis();
        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        int nGroupsToFilterGenotypingDataOn = VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gsver, true).size();

        Map<String, Collection<String>> individualsByPop = new HashMap<>();
        Map<String, HashMap<String, Float>> annotationFieldThresholdsByPop = new HashMap<>();
        List<List<String>> callsetIds = gsver.getAllCallSetIds();
        for (int i = 0; i < callsetIds.size(); i++) {
            individualsByPop.put(gsver.getGroupName(i), callsetIds.get(i).isEmpty() ? MgdbDao.getProjectIndividuals(sModule, projId) /* no selection means all selected */ : callsetIds.get(i).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet()));
            annotationFieldThresholdsByPop.put(gsver.getGroupName(i), gsver.getAnnotationFieldThresholds(i));
        }

        Collection<String> individualsToExport = gsver.getExportedIndividuals().size() > 0 ? gsver.getExportedIndividuals() : MgdbDao.getProjectIndividuals(sModule, projId);

        long count = countVariants(gsver, true);
        MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(sModule, token, false, false, false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
        if (nGroupsToFilterGenotypingDataOn > 0 && nTempVarCount == 0)
        {
            progress.setError(MgdbDao.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return;
        }

        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gsver/*, getSequenceIDsBeingFilteredOn(gsver.getRequest().getSession(), sModule)*/, true);
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        final BasicDBList variantQueryDBList = variantDataQueries.size() == 1 ? variantDataQueries.iterator().next() : new BasicDBList();

        Document variantQuery = nTempVarCount == 0 && !variantDataQueries.isEmpty() ? new Document("$and", variantQueryDBList) : new Document();
        String usedVarCollName = nTempVarCount == 0 ? mongoTemplate.getCollectionName(VariantData.class) : tmpVarColl.getNamespace().getCollectionName();

        Authentication auth = tokenManager.getAuthenticationFromToken(tokenManager.readToken(gsver.getRequest()));
        if (gsver.shallApplyMatrixSizeLimit())
        {    // make sure the matrix is not too big
            int nMaxBillionGenotypesInvolved = 1;    // default
            try {
                nMaxBillionGenotypesInvolved = Integer.parseInt(appConfig.get("maxExportableBillionGenotypes_" + (auth == null ? "anonymousUser" : auth.getName())));
            }
            catch (Exception ignored1) {
                   try {
                    nMaxBillionGenotypesInvolved = Integer.parseInt(appConfig.get("maxExportableBillionGenotypes"));
                }
                catch (Exception ignored2)
                {}
            }

            if (nMaxBillionGenotypesInvolved == 0)
            {
                progress.setError("You are not allowed to export any genotyping data.");
                return;
            }

            BigInteger matrixSize = BigInteger.valueOf(mongoTemplate.getCollection(usedVarCollName).countDocuments(variantQuery)).multiply(BigInteger.valueOf(individualsToExport.size()));
            BigInteger maxAllowedSize = BigInteger.valueOf(1000000000).multiply(BigInteger.valueOf(nMaxBillionGenotypesInvolved));

            if (matrixSize.divide(maxAllowedSize).intValue() >= 1)
            {
                progress.setError("You may only export up to " + nMaxBillionGenotypesInvolved + " billion genotypes. The current selection contains " + BigDecimal.valueOf(matrixSize.longValue()).divide(BigDecimal.valueOf(1000000000)).setScale(2, BigDecimal.ROUND_HALF_UP) + " billion genotypes.");
                return;
            }
        }

        progress.addStep("Identifying matching variants");

        OutputStream os = null;

        try
        {
        	String username = AbstractTokenManager.getUserNameFromAuthentication(tokenManager.getAuthenticationFromToken(token));

            AbstractIndividualOrientedExportHandler individualOrientedExportHandler = AbstractIndividualOrientedExportHandler.getIndividualOrientedExportHandlers().get(gsver.getExportFormat());
            AbstractMarkerOrientedExportHandler markerOrientedExportHandler = AbstractMarkerOrientedExportHandler.getMarkerOrientedExportHandlers().get(gsver.getExportFormat());

            String filename = sModule + "__project" + projId + "__" + new SimpleDateFormat("yyyy-MM-dd").format(new Date()) + "__" + count + "variants__" + gsver.getExportFormat().replace(".", "_") + "." + (individualOrientedExportHandler != null ? individualOrientedExportHandler : markerOrientedExportHandler).getExportArchiveExtension();

            LOG.info((gsver.isKeepExportOnServer() ? "On-server" : "Direct-download") + " export requested: " + processId);

            String relativeOutputFolder = FRONTEND_URL + File.separator + (!gsver.isKeepExportOnServer() ? TMP_OUTPUT_DDL_FOLDER : TMP_OUTPUT_FOLDER) + File.separator + username + File.separator + Helper.convertToMD5(processId) + File.separator;
            File outputLocation = new File(servletContext.getRealPath(File.separator + relativeOutputFolder));
            if (!outputLocation.exists() && !outputLocation.mkdirs()) {
                throw new Exception("Unable to create folder: " + outputLocation);
            }
            
            HttpSession session = gsver.getRequest().getSession();	// get it before starting to use the response (otherwise CORS invocations may make it invalid when we try and use it) 
            os = new FileOutputStream(new File(outputLocation.getAbsolutePath() + File.separator + filename));
            response.setContentType("text/plain");
            String exportURL = gsver.getRequest().getContextPath() + "/" + relativeOutputFolder.replace(File.separator, "/") + filename;
            LOG.debug("Export file for process " + processId + ": " + exportURL);
            response.getWriter().write(exportURL);
            response.flushBuffer();

            GenotypingProject project = mongoTemplate.findById(projId, GenotypingProject.class);
            Map<String, InputStream> readyToExportFiles = new HashMap<>();
            String sCitingText = appConfig.get("howToCite");
            if (sCitingText == null)
                sCitingText = "Please cite Gigwa as follows:\nGuilhem Sempéré, Adrien Pétel, Mathieu Rouard, Julien Frouin, Yann Hueber, Fabien De Bellis, Pierre Larmande,\nGigwa v2—Extended and improved genotype investigator, GigaScience, Volume 8, Issue 5, May 2019, giz051, https://doi.org/10.1093/gigascience/giz051";
            String projDesc = project.getDescription();
            if (projDesc != null && projDesc.contains("HOW TO CITE"))
                sCitingText += (sCitingText.length() > 0 ? "\n\n" : "") + "Please cite project data as follows:\n" + projDesc.substring(projDesc.indexOf("HOW TO CITE") + 11).replaceAll("\n\n*", "\n").trim();
            if (sCitingText.length() > 0)
                readyToExportFiles.put("HOW_TO_CITE.txt", new ByteArrayInputStream(sCitingText.getBytes("UTF-8")));

            final OutputStream finalOS = os;
            ArrayList<GenotypingSample> samplesToExport = MgdbDao.getSamplesForProject(sModule, projId, individualsToExport);
            final Integer nAssembly = Assembly.getThreadBoundAssembly();
            if (individualOrientedExportHandler != null)
            {
                if (!progress.isAborted()) {
                    Thread exportThread = new SessionAttributeAwareThread(session) {
                        public void run() {
                        	Assembly.setThreadAssembly(nAssembly);	// set it once and for all
                            try {
                                ExportOutputs exportOutputs = individualOrientedExportHandler.createExportFiles(sModule, Assembly.getThreadBoundAssembly(), nTempVarCount == 0 ? null : usedVarCollName, nTempVarCount == 0 ? variantQueryDBList : new BasicDBList(), count, processId, individualsByPop, annotationFieldThresholdsByPop, samplesToExport, progress);

                                for (String step : individualOrientedExportHandler.getStepList())
                                    progress.addStep(step);
                                progress.moveToNextStep();
                                individualOrientedExportHandler.exportData(finalOS, sModule, Assembly.getThreadBoundAssembly(), AbstractTokenManager.getUserNameFromAuthentication(auth), exportOutputs, true, progress, nTempVarCount == 0 ? null : usedVarCollName, varQueryWrapper, count, null, gsver.getMetadataFields(), IExportHandler.getIndividualPopulations(individualsByPop, true), readyToExportFiles);
                                if (!progress.isAborted()) {
                                    LOG.info("exportVariants (" + gsver.getExportFormat() + ") took " + (System.currentTimeMillis() - before) / 1000d + "s to process " + count + " variants and " + individualsToExport.size() + " individuals");
                                    progress.markAsComplete();
                                }
                            }
                            catch (Exception e) {
                                LOG.error("Error exporting data", e);
                                progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
                            }
                            finally {
                                try
                                {
                                    finalOS.close();
                                }
                                catch (IOException ignored)
                                {}
                            }
                        }
                    };
                    exportThread.start();
                }
            }
            else if (markerOrientedExportHandler != null)
            {
                for (String step : markerOrientedExportHandler.getStepList()) {
                    progress.addStep(step);
                }
                progress.moveToNextStep();    // done with identifying variants

                String contentType = markerOrientedExportHandler.getExportContentType();
                if (contentType != null && contentType.trim().length() > 0)
                    response.setContentType(contentType);

                Thread exportThread = new SessionAttributeAwareThread(session) {
                    public void run() {
                    	Assembly.setThreadAssembly(nAssembly);	// set it once and for all
                        try {
                            markerOrientedExportHandler.exportData(finalOS, sModule, Assembly.getThreadBoundAssembly(), AbstractTokenManager.getUserNameFromAuthentication(auth), progress, nTempVarCount == 0 ? null : usedVarCollName, varQueryWrapper, count, null, individualsByPop, annotationFieldThresholdsByPop, samplesToExport, gsver.getMetadataFields(), null);
                            if (!progress.isAborted() && progress.getError() == null) {
                                LOG.info("exportVariants (" + gsver.getExportFormat() + ") took " + (System.currentTimeMillis() - before) / 1000d + "s to process " + count + " variants and " + individualsToExport.size() + " individuals");
                                progress.markAsComplete();
                            }
                        }
                        catch (Exception e) {
                            LOG.error("Error exporting data", e);
                            progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
                        }
                        finally {
                            try
                            {
                                finalOS.close();
                            }
                            catch (IOException ignored)
                            {}
                        }
                    }
                };
                exportThread.start();
            }
            else
                throw new Exception("No export handler found for format " + gsver.getExportFormat());
        } catch (Throwable t) {
            LOG.error("Error exporting data", t);
            progress.setError("Error exporting data: " + t.getClass().getSimpleName() + (t.getMessage() != null ? " - " + t.getMessage() : ""));
            return;
        }
    }

//    @Override
//    public int getSequenceFilterCount(HttpServletRequest request, String sModule) throws IOException {
//        int result = -1;
//        File sequenceListFile = new File(servletContext.getRealPath(SEQLIST_FOLDER + File.separator + request.getSession().getId() + "_" + sModule));
//        if (sequenceListFile.exists() && sequenceListFile.length() > 0) {
//            try (LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(sequenceListFile))) {
//                lineNumberReader.skip(Long.MAX_VALUE);
//                result = lineNumberReader.getLineNumber();
//            }
//        }
//        return result;
//    }
//
//    @Override
//    public ArrayList<String> getSequenceIDsBeingFilteredOn(HttpSession session, String sModule) {
//        ArrayList<String> sequences = new ArrayList<>();
//        File selectionFile = new File(session.getServletContext().getRealPath(SEQLIST_FOLDER) + File.separator + session.getId() + "_" + sModule);
//        if (selectionFile.exists() && selectionFile.length() > 0) {
//            Scanner sc = null;
//            try {
//                sc = new Scanner(selectionFile);
//                sc.nextLine();    // skip queryKey line
//                while (sc.hasNextLine()) {
//                    sequences.add(sc.nextLine().trim());
//
//                } // skip queryKey line
//            } catch (FileNotFoundException ex) {
//                LOG.debug("couldn't find sequence list file", ex);
//            } finally {
//                sc.close();
//            }
//        }
//        return sequences.isEmpty() ? null : sequences;
//    }
//
//    @Override
//    public void clearSequenceFilterFile(HttpServletRequest request, String sModule) {
//        File selectionFile = new File(servletContext.getRealPath(SEQLIST_FOLDER) + File.separator + request.getSession().getId() + "_" + sModule);
//        if (selectionFile.exists()) {
//            selectionFile.delete();
//        }
//    }
//
//    @Override
//    public String getSequenceFilterQueryKey(HttpServletRequest request, String sModule) throws IOException {
//
//        String qk = null;
//        File sequenceListFile = new File(servletContext.getRealPath(SEQLIST_FOLDER + File.separator + request.getSession().getId() + "_" + sModule));
//        if (sequenceListFile.exists() && sequenceListFile.length() > 0) {
//            try (LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(sequenceListFile))) {
//                qk = lineNumberReader.readLine();
//            }
//        }
//        return qk;
//    }
    
    static public void cleanupExpiredExportData(ServletContext sc) throws IOException {
        Map<String, Long> folderToDelayMap = new LinkedHashMap<>() {{
        	put(TMP_OUTPUT_DDL_FOLDER, DDL_EXPORT_EXPIRATION_DELAY_MILLIS);
        	put(TMP_OUTPUT_FOLDER, EXPORT_EXPIRATION_DELAY_MILLIS);
        	put(TMP_OUTPUT_EXTRACTION_FOLDER, DDL_EXPORT_EXPIRATION_DELAY_MILLIS);
        }};
        
        long nowMillis = new Date().getTime();
        for (Entry<String, Long> entry : folderToDelayMap.entrySet()) {
            File filterOutputLocation = new File(sc.getRealPath(FRONTEND_URL + File.separator + entry.getKey()));
            if (filterOutputLocation.exists() && filterOutputLocation.isDirectory()) {
            	for (File userFolder : filterOutputLocation.listFiles())
            		if (userFolder.isDirectory()) {
                    	for (File processFolder : userFolder.listFiles())
                    		if (processFolder.isDirectory()) {
                    			if (nowMillis - processFolder.lastModified() > entry.getValue()) {
                					FileSystemUtils.deleteRecursively(processFolder);
//                					System.err.println("deleting expired folder: " + processFolder.getAbsolutePath());
                    			}
                    			else {
	                    			for (File exportFileOrFolder : processFolder.listFiles()) {
	                    				if (nowMillis - exportFileOrFolder.lastModified() > entry.getValue()) {
		                    				if (exportFileOrFolder.isDirectory()) {
		                    					FileSystemUtils.deleteRecursively(exportFileOrFolder);
//		                    					System.err.println("deleting expired folder: " + exportFileOrFolder.getAbsolutePath());
		                    				}
		                    				else {
		                    					exportFileOrFolder.delete();
//		                    					System.err.println("deleting expired file: " + exportFileOrFolder.getAbsolutePath());
		                    				}
	                    				}
	                    			}
	                    			if (processFolder.listFiles().length == 0) {
	                					FileSystemUtils.deleteRecursively(processFolder);
//	                					System.err.println("deleting empty folder: " + processFolder.getAbsolutePath());
	                    			}
                    			}
                    		}

            			if (userFolder.listFiles().length == 0 && !TMP_OUTPUT_EXTRACTION_FOLDER.equals(entry.getKey()) /* avoid interfering with ongoing exports */) {
        					FileSystemUtils.deleteRecursively(userFolder);
//        					System.err.println("deleting empty folder: " + userFolder.getAbsolutePath());
            			}
            		}
            }
        }
    }

    @Override
    public boolean abortProcess(String processID) {
        ProgressIndicator progress = ProgressIndicator.get(processID);
        if (progress != null) {
            progress.abort();
            LOG.debug("Aborting process: " + processID + " [" + progress.hashCode() + "]");
            return true;
        }
        return false;
    }

    @Override
    public void dropTempColl(String sModule, String processID) throws InterruptedException {
        String collName = MongoTemplateManager.getTemporaryVariantCollection(sModule, processID, false, false, false).getNamespace().getCollectionName();
        MongoTemplateManager.get(sModule).dropCollection(collName);
        LOG.debug("Dropped collection " + sModule + "." + collName);
    }

    @Override
    public Collection<String> distinctSequencesInSelection(HttpServletRequest request, String sModule, int projId, String processID) throws InterruptedException {
        String sShortProcessID = processID/*.substring(1 + processID.indexOf('|'))*/;
        MongoCollection<Document> tmpVarColl = MongoTemplateManager.getTemporaryVariantCollection(sModule, sShortProcessID, false, false, false);
        if (tmpVarColl.estimatedDocumentCount() == 0) {
            return listSequences(request, sModule, projId);    // working on full dataset
        }
        List<String> distinctSequences = tmpVarColl.distinct(Assembly.getThreadBoundVariantRefPosPath() + "." + ReferencePosition.FIELDNAME_SEQUENCE, String.class).into(new ArrayList<>());
        TreeSet<String> sortedResult = new TreeSet<>(new AlphaNumericComparator());
        sortedResult.addAll(distinctSequences);
        return sortedResult;
    }


    @Override
    public String getQueryKey(MgdbSearchVariantsRequest gsvr) {
        String info[] = Helper.getInfoFromId(gsvr.getVariantSetId(), 2);
        int projId = Integer.parseInt(info[1]);
        String queryKey = projId + ":" + Assembly.getThreadBoundAssembly() + ":"
                        + gsvr.getSelectedVariantTypes() + ":"
                        + gsvr.getReferenceName() + ":"
                        + (gsvr.getStart() == null ? "" : gsvr.getStart()) + ":"
                        + (gsvr.getEnd() == null ? "" : gsvr.getEnd()) + ":"
                        + gsvr.getAlleleCount() + ":"
                        + gsvr.getGeneName() + ":"
                        + gsvr.getSelectedVariantIds() + ":";
        
        List<List<String>> callsetIds = gsvr.getAllCallSetIds();
        for (int i = 0; i < callsetIds.size(); i++) {
            queryKey += callsetIds.get(i) + ":"
                        + gsvr.getAnnotationFieldThresholds(i) + ":"
                        + gsvr.getGtPattern(i) + ":"
                        + gsvr.getMostSameRatio(i) + ":"
                        + gsvr.getMinMissingData(i) + ":"
                        + gsvr.getMaxMissingData(i) + ":"
                        + gsvr.getMinHeZ(i) + ":"
                        + gsvr.getMaxHeZ(i) + ":"
                        + gsvr.getMinMaf(i) + ":"
                        + gsvr.getMaxMaf(i) + ":";
        }
        queryKey += gsvr.getDiscriminate() + ":"
                  + gsvr.getVariantEffect();
//        System.out.println(queryKey);
        return Helper.convertToMD5(queryKey);
    }

    /**
     * Import a sequence from a fasta in the database
     *
     * @param module
     * @param filePath fasta file local or remote path (accept ftp:// )
     * @param mode
     * @return boolean true if import succeded
     * @throws java.lang.Exception
     */
    public boolean importSequenceInDB(String module, String filePath, String mode) throws Exception {
        boolean success = false;
        try {
            SequenceImport.main(new String[]{module, filePath, mode});
            success = true;
        } catch (Exception ex) {
            throw ex;
        }
        return success;
    }

    /**
     * get all run in a project
     *
     * @param module
     * @param projId
     * @return
     */
    public List<String> getRunList(String module, int projId) {

        List<String> listRun;
        MongoTemplate mongoTemplate = MongoTemplateManager.get(module);
        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_RUNS);
           q.addCriteria(Criteria.where("_id").is(projId));
        GenotypingProject project = mongoTemplate.findOne(q, GenotypingProject.class);
        listRun = project.getRuns();
        if (listRun == null) {
            return new ArrayList<>();
        } else {
            return listRun;
        }
    }

    /**
     * get information on a list of sequence
     *
     * @param module
     * @param sequenceList
     * @return
     */
    public Map<String, Map<String, Object>> getSequenceInfo(String module, List<String> sequenceList) {

        Map<String, Map<String, Object>> listSeqInfo = new HashMap<>();
        Map<String, Object> info = new HashMap<>();
        Document sequence;

        // no need to filter on projId since we are searching on sequenceId
        ArrayList<BasicDBObject> aggregationParam = new ArrayList<>();
        aggregationParam.add(new BasicDBObject("$match", new BasicDBObject("_id", new BasicDBObject("$in", sequenceList))));

        MongoCursor<Document> genotypingDataCursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(aggregationParam).allowDiskUse(true).iterator();

        if (genotypingDataCursor != null && genotypingDataCursor.hasNext()) {

            while (genotypingDataCursor.hasNext()) {

                sequence = genotypingDataCursor.next();

                // empty the table
                info.clear();

                // no need for sequence?
                info.put("sequence", sequence.get(Sequence.FIELDNAME_SEQUENCE));

                info.put("length", sequence.get(Sequence.FIELDNAME_LENGTH));
                info.put("checksum", sequence.get(Sequence.FIELDNAME_CHECKSUM));

                listSeqInfo.put((String) sequence.get("_id"), info);

            }
        }
        return listSeqInfo;
    }

    /**
     * get a list of variant in ga4gh format from a MongoCursor<Document>
     *
     * @param module
     * @param projId
     * @param cursor
     * @param samples
     * @return List<Variant>
     * @throws AvroRemoteException
     */
    public List<Variant> getVariantListFromDBCursor(String module, int projId, MongoCursor<Document> cursor, Collection<GenotypingSample> samples)
    {
//        long before = System.currentTimeMillis();
        LinkedHashMap<Comparable, Variant> varMap = new LinkedHashMap<>();

        String refPosPath = Assembly.getThreadBoundVariantRefPosPath();
        
        // parse the cursor to create all GAVariant
        while (cursor.hasNext()) {
            Document obj = cursor.next();
            // save the Id of each variant in the cursor
            String id = (String) obj.get("_id");
            List<String> knownAlleles = ((List<String>) obj.get(VariantData.FIELDNAME_KNOWN_ALLELES));

            Variant.Builder variantBuilder = Variant.newBuilder().setId(Helper.createId(module, projId, id.toString())).setVariantSetId(Helper.createId(module, projId));

            Document rp = (Document) Helper.readPossiblyNestedField(obj, refPosPath, "; ", null);
            if (rp == null)
                variantBuilder.setReferenceName("").setStart(0).setEnd(0);
            else {
                String chr = (String) rp.get(ReferencePosition.FIELDNAME_SEQUENCE);
                variantBuilder.setReferenceName(chr != null ? chr : "");

                Long start = (Long) rp.get(ReferencePosition.FIELDNAME_START_SITE);
                variantBuilder.setStart(start != null ? start : 0);

                Long end = (Long) rp.get(ReferencePosition.FIELDNAME_END_SITE);
                if (end == null && start != null)
                	end = start;
                variantBuilder.setEnd(end != null ? end : 0);
            }
            if (knownAlleles != null && knownAlleles.size() > 0) {
	            variantBuilder.setReferenceBases(knownAlleles.get(0)); // reference is the first one in VCF files
	            variantBuilder.setAlternateBases(knownAlleles.subList(1, knownAlleles.size()));
            }
            else
                variantBuilder.setReferenceBases("");	// the DTO used here requires a reference allele to be provided: the chosen convention is to provide a single allele coded as an empty string when no allele is known for a variant

            // add the annotation map to the variant
            Map<String, List<String>> annotations = new HashMap<>();
            List<String> infoType = new ArrayList<>();
            infoType.add((String) obj.get(VariantData.FIELDNAME_TYPE));
            annotations.put("type", infoType);
            variantBuilder.setInfo(annotations);
            
            varMap.put(id, variantBuilder.build());
        }

        // get the VariantRunData containing annotations
        ArrayList<BasicDBObject> pipeline = new ArrayList<>();
        // wanted fields
        BasicDBObject fields = new BasicDBObject();
        fields.put(VariantRunData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_GENE, 1);
        fields.put(VariantRunData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME, 1);

        // get the genotype for wanted individuals/callSet only
        boolean fGotMultiSampleIndividuals = false;
        HashSet<String> involvedIndividuals = new HashSet<>();
        for (GenotypingSample sample : samples){
            fields.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + sample.getId(), 1);
            if (!involvedIndividuals.add(sample.getIndividual()))
            	fGotMultiSampleIndividuals = true;
        }

        BasicDBList matchAndList = new BasicDBList();
        matchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", varMap.keySet())));
        matchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID, projId));
        if (!samples.isEmpty())
            matchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_RUNNAME, new BasicDBObject("$in", samples.stream().map(sp -> sp.getRun()).distinct().collect(Collectors.toList()))));
        pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", matchAndList)));
        pipeline.add(new BasicDBObject("$project", fields));
        if (samples.isEmpty())    // if no genotypes are expected back then we assume we're building the result table (thus we need to include variant name & effect when available in one of then runs)
            pipeline.add(new BasicDBObject("$sort", new BasicDBObject(AbstractVariantData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME, -1)));  // if some VariantRunData records have gene info they will appear first, which will make that info available for building the result table

        HashSet<String> variantsForWhichAnnotationWasRetrieved = new HashSet<>();

        MongoCursor<Document> genotypingDataCursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(pipeline).allowDiskUse(true).iterator();
        if (!genotypingDataCursor.hasNext())
        	for (Comparable varId : varMap.keySet()) {	// create empty Call documents for all requested variants
                Variant var = varMap.get(varId);
                TreeSet<Call> calls = new TreeSet(new AlphaNumericComparator<Call>());    // for automatic sorting
                Builder emptyCall = Call.newBuilder().setGenotype(new ArrayList<>());
        		for (GenotypingSample sample : samples) {
        			emptyCall.setCallSetId(Helper.createId(module, projId, sample.getIndividual()));
                    calls.add(emptyCall.build());
	        	}
            	var.setCalls(new ArrayList<Call>(calls));	// add the call list
	        }
        else while (genotypingDataCursor.hasNext()) {
            Document variantObj = genotypingDataCursor.next();
            String varId = (String) Helper.readPossiblyNestedField(variantObj, "_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, "; ", null);
            Variant var = varMap.get(varId);
            if (var == null /* should not happen! */|| variantsForWhichAnnotationWasRetrieved.contains(varId))
                continue;

            variantsForWhichAnnotationWasRetrieved.add(varId);
            TreeSet<Call> calls = new TreeSet(new AlphaNumericComparator<Call>());    // for automatic sorting

            // for each annotation field
            for (String key : variantObj.keySet()) {
                switch (key) {
                    // this goes in Call  || should not be called if sp field is not present
                    case VariantRunData.FIELDNAME_SAMPLEGENOTYPES:
                        // get genotype map
                        Map<String, Object> callMap = (Map<String, Object>) variantObj.get(key);

                        // for each individual/CallSet
                        for (GenotypingSample sample : samples) {
                            Document callObj = (Document) callMap.get("" + sample.getId());
                            double[] gl;
                            List<Double> listGL = new ArrayList<>();

                            List<Integer> genotype = new ArrayList<>();
                            Map<String, List<String>> aiCall = new HashMap<>();
                            String phaseSet = null;

                            if (callObj != null)
                            {
                                Map<String, Object> callAdditionalInfo = (Map<String, Object>) callObj.get("ai");

                                // if field ai is present
                                if (callAdditionalInfo != null)
                                    for (String aiKey : callAdditionalInfo.keySet()) {
                                        if (aiKey.equals(VCFConstants.GENOTYPE_PL_KEY))
                                        {
                                            gl = GenotypeLikelihoods.fromPLField(callAdditionalInfo.get(aiKey).toString()).getAsVector();
                                            for (int h = 0; h < gl.length; h++)
                                                listGL.add(gl[h]);
                                        }
                                        switch (aiKey) {
                                            case VCFConstants.GENOTYPE_LIKELIHOODS_KEY:
                                                listGL = (List<Double>) callAdditionalInfo.get(aiKey);
                                                break;

                                            case VCFConstants.PHASE_SET_KEY:
                                            case VariantData.GT_FIELD_PHASED_ID:
                                                phaseSet = callAdditionalInfo.get(aiKey).toString();
                                                break;

                                            default:
                                            	Object val = callAdditionalInfo.get(aiKey);
                                            	if (val != null)
	                                                aiCall.put(aiKey, Arrays.asList(val.toString()));
                                                break;
                                        }
                                    }

                                // get GT info
                                String gt = (String) callObj.get("gt");

                                if (gt == null || gt.startsWith(".")) {
                                    // if we don't know the genotype, do nothing
                                } else {
                                    String[] gen;
                                    if (gt.contains("/")) {
                                        gen = gt.split("/");
                                    } else {
                                        gen = gt.split(Helper.ID_SEPARATOR);
                                    }
                                    for (String gen1 : gen) {
                                        genotype.add(Integer.parseInt(gen1));
                                    }
                                }
                            }
                            
                            if (fGotMultiSampleIndividuals)
                            	aiCall.put("sample", Arrays.asList("" + sample.getSampleName()));
                            Call call = Call.newBuilder()
                                    .setCallSetId(Helper.createId(module, projId, sample.getIndividual()))
                                    .setGenotype(genotype)
                                    .setGenotypeLikelihood(listGL)
                                    .setPhaseset(phaseSet)
                                    .setInfo(aiCall)
                                    .build();

                            calls.add(call);
                        }
                        break;

                    case VariantRunData.SECTION_ADDITIONAL_INFO:
                        Map<String, Object> additionalInfos = (Map<String, Object>) variantObj.get(key);
                        for (String subKey : additionalInfos.keySet()) {

                            if (subKey.equals("") || subKey.equals(VcfImport.ANNOTATION_FIELDNAME_ANN) || subKey.equals(VcfImport.ANNOTATION_FIELDNAME_CSQ)) {
                                // if VCF has empty field (";") do not retrieve it

                                // field EFF should be stored in variantAnnotation !
                                // stored in ai for the moment, not supported by ga4gh
                                // ANN (vcf 4.2) is stored in variantAnnotation
                            } else if (subKey.equals(VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_GENE)) {
                                List<String> listGene = (List<String>) additionalInfos.get(subKey);
                                var.getInfo().put(subKey, listGene);
                            } else if (subKey.equals(VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME)) {
                                List<String> listEffect = (List<String>) additionalInfos.get(subKey);
                                var.getInfo().put(subKey, listEffect);
                            } else {

                            }
                        }
                        break;
                    default:
                        // "_id" and "_class", do nothing
                        break;
                }
            }
            var.setCalls(new ArrayList<Call>(calls));	// add the call list
        }

//        LOG.debug("getVariantListFromDBCursor took " + (System.currentTimeMillis() - before) / 1000f + "s for " + varMap.size() + " variants and " + samples.size() + " samples");
        return new ArrayList<Variant>(varMap.values());
    }

    /**
     * get the Metadata number from a VCFCompoundLine
     *
     * @param vcf
     * @return String number
     */
    private String getNumber(VCFCompoundHeaderLine vcf) {

        String number;
        switch (vcf.getCountType()) {
            case A:
                number = "";
                break;
            case G:
                number = "";
                break;
            case INTEGER:
                number = Integer.toString(vcf.getCount());
                break;
            case R:
                number = "";
                break;
            case UNBOUNDED:

                // ga4gh python serveur return "." when unbounded
                // but no information about it in ga4gh documentation
                number = Integer.toString(-1);
                break;
            default:
                number = "";
                break;
        }
        return number;
    }

    /**
     * return the progress
     *
     * @param processID
     * @return ProgressIndicator
     */
    public ProgressIndicator getProgressIndicator(String processID) {
        return ProgressIndicator.get(processID/*.substring(1 + processID.indexOf('|'))*/);
    }

    /**
     * get the metadata for a Project/VariantSet
     *
     * @param module
     * @param proj
     * @return List<List<VariantSetMetadata>>
     */
    private List<VariantSetMetadata> getMetadataList(String module, String proj) {

        List<VariantSetMetadata> listMetadata;
        VariantSetMetadata metadata;
        Map<String, VariantSetMetadata> mapMetadata = new HashMap<>();
        DBVCFHeader vcfHeader;

        VCFInfoHeaderLine vci;
        VCFFilterHeaderLine vcf;
        VCFFormatHeaderLine vcfo;
        VCFSimpleHeaderLine vcm;
        VCFHeaderLine vcmo;
        String metadataKey;
        int i;
        Map<String, List<String>> info = new HashMap<>();

        // get the list of vcf header for this project
        BasicDBObject whereQuery = new BasicDBObject();
        whereQuery.put("_id." + DBVCFHeader.VcfHeaderId.FIELDNAME_PROJECT, Integer.parseInt(proj));
        MongoCursor<Document> cursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class)).find(whereQuery).iterator();

        if (cursor != null && cursor.hasNext()) {

            while (cursor.hasNext()) {

                // get the vcf header of each run of a project
                vcfHeader = DBVCFHeader.fromDocument(cursor.next());

                // used to create id for each row header
                i = 0;

                // fill metadata list
                for (String key : vcfHeader.getmInfoMetaData().keySet()) {

                    vci = vcfHeader.getmInfoMetaData().get(key);

                    // store the key to make sure a header is only present once
                    // (if project has multiple run)
                    metadataKey = vci.getKey() + "." + vci.getID();

                    metadata = VariantSetMetadata.newBuilder()
                            .setDescription(vci.getDescription())
                            .setType(vci.getType().toString())
                            .setValue(vci.getValue())
                            .setId(Integer.toString(i))
                            .setKey(metadataKey)
                            .setNumber(getNumber(vci))
                            .setInfo(info)
                            .build();
                    i++;
                    mapMetadata.put(metadataKey, metadata);
                }
                for (String key : vcfHeader.getmFilterMetaData().keySet()) {

                    vcf = vcfHeader.getmFilterMetaData().get(key);
                    metadataKey = vcf.getKey() + "." + vcf.getID();

                    metadata = VariantSetMetadata.newBuilder()
                            .setId(Integer.toString(i))
                            .setType("String")
                            .setKey(metadataKey)
                            .setValue(vcf.getValue())
                            .setInfo(info)
                            .setDescription("")
                            .setNumber("")
                            .build();
                    i++;
                    mapMetadata.put(metadataKey, metadata);

                }
                for (String key : vcfHeader.getmFormatMetaData().keySet()) {

                    vcfo = vcfHeader.getmFormatMetaData().get(key);

                    metadataKey = vcfo.getKey() + "." + vcfo.getID();

                    metadata = VariantSetMetadata.newBuilder()
                            .setId(Integer.toString(i))
                            .setKey(metadataKey)
                            .setValue(vcfo.getValue())
                            .setType(vcfo.getType().toString())
                            .setDescription(vcfo.getDescription())
                            .setInfo(info)
                            .setNumber(getNumber(vcfo))
                            .build();
                    i++;
                    mapMetadata.put(metadataKey, metadata);

                }
                for (String key : vcfHeader.getmMetaData().keySet()) {

                    vcm = vcfHeader.getmMetaData().get(key);
                    metadataKey = vcm.getKey() + "." + vcm.getID();

                    metadata = VariantSetMetadata.newBuilder()
                            .setId(Integer.toString(i))
                            .setKey(metadataKey)
                            .setValue(vcm.getValue())
                            .setDescription("")
                            .setNumber("")
                            .setType("String")
                            .setInfo(info)
                            .build();
                    i++;
                    mapMetadata.put(metadataKey, metadata);

                }

                for (String key : vcfHeader.getmOtherMetaData().keySet()) {

                    vcmo = vcfHeader.getmOtherMetaData().get(key);
                    metadataKey = vcmo.getKey();

                    metadata = VariantSetMetadata.newBuilder()
                            .setId(Integer.toString(i))
                            .setKey(metadataKey)
                            .setValue(vcmo.getValue())
                            .setDescription("")
                            .setType("String")
                            .setNumber("")
                            .setInfo(info)
                            .build();
                    i++;
                    mapMetadata.put(metadataKey, metadata);
                }
            }
            cursor.close();
        }
        // get a list of variantSetMetadata from the map
        listMetadata = new ArrayList<>(mapMetadata.values());

        return listMetadata;
    }

    /*
     * GA4GH methods, from methods interface v0.6.1 - opencb/ga4gh 04/2016
     */
    @Override
    public VariantSet getVariantSet(String id) throws AvroRemoteException
    {
        VariantSet variantSet = null;
        // get information from id
        String[] info = Helper.getInfoFromId(id, 2);
        if (info != null)
            try
            {
                String module = info[0];
                String projId = info[1];

                Query q = new Query(Criteria.where("_id").is(Integer.parseInt(projId)));
                q.fields().include(GenotypingProject.FIELDNAME_NAME);

                GenotypingProject proj = MongoTemplateManager.get(module).findOne(q, GenotypingProject.class);
                List<VariantSetMetadata> metadata = getMetadataList(module, projId);
                if (proj.getDescription() != null)
                {
                    VariantSetMetadata vsmd = new VariantSetMetadata();
                    vsmd.setKey(AbstractVariantData.VCF_CONSTANT_DESCRIPTION);
                    vsmd.setValue(proj.getDescription());
                    metadata.add(vsmd);
                }
                // Checksum : MD5 of the upper-case sequence excluding all whitespace characters
                // length == 0 since we don't have this information in VCF files
                variantSet = VariantSet.newBuilder()
                        .setId(id)
                        .setDatasetId(module)
                        .setName(proj.getName())
                        .setReferenceSetId(module)
                        .setMetadata(metadata)
                        .build();
            }
            catch (NumberFormatException nfe)
            {}
        return variantSet;
    }

    @Override
    public Variant getVariant(String id) throws AvroRemoteException {
        return getVariantWithGenotypes(id, new ArrayList<>() /* all individuals */);
    }

    public Variant getVariantWithGenotypes(String id, Collection<String> listInd) throws NumberFormatException, AvroRemoteException {
        String[] info = id.split(Helper.ID_SEPARATOR);
        MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
        MongoCursor<Document> cursor = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantData.class)).find(new BasicDBObject("_id", info[2])).iterator();

        Variant variant = null;
        if (cursor != null && cursor.hasNext()) {
            List<Criteria> sampleQueryCriteria = new ArrayList<>();
            sampleQueryCriteria.add(Criteria.where(GenotypingSample.FIELDNAME_PROJECT_ID).is(Integer.parseInt(info[1])));
            if (!listInd.isEmpty())
                sampleQueryCriteria.add(Criteria.where(GenotypingSample.FIELDNAME_INDIVIDUAL).in(listInd));
            if (info.length == 4)	// run id may optionally be appended to variant id, to restrict samples to those involved in the run
                sampleQueryCriteria.add(Criteria.where(GenotypingSample.FIELDNAME_RUN).is(info[3]));
            Collection<GenotypingSample> samples = mongoTemplate.find(new Query(new Criteria().andOperator(sampleQueryCriteria.toArray(new Criteria[sampleQueryCriteria.size()]))), GenotypingSample.class);
            variant = getVariantListFromDBCursor(info[0], Integer.parseInt(info[1]), cursor, samples).get(0);
            cursor.close();
        }
        return variant;
    }

    @Override
    public CallSet getCallSet(String id) throws AvroRemoteException {
        CallSet callSet = null;

        // get information from id
        String[] info = Helper.getInfoFromId(id, 3);
        if (info == null) {

            // wrong number of param or wrong module name
        } else {
            String module = info[0];
            int projId = Integer.parseInt(info[1]);
            String name = info[2];

            List<String> listVariantSetId = new ArrayList<>();
            listVariantSetId.add(Helper.createId(module, info[1]));

            try {
                // check if the callSet is in the list
                if (MgdbDao.getProjectIndividuals(module, projId).contains(name))
                    callSet = CallSet.newBuilder().setId(id).setName(name).setVariantSetIds(listVariantSetId).setSampleId(null).build();
            } catch (ObjectNotFoundException ex) {
                java.util.logging.Logger.getLogger(GigwaGa4ghServiceImpl.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return callSet;
    }

    @Override
    public ReferenceSet getReferenceSet(String id) throws AvroRemoteException {
        ReferenceSet referenceSet = null;

        MongoTemplate mongoTemplate = MongoTemplateManager.get(id);
        if (mongoTemplate == null) {

        } else {
            List<String> list = new ArrayList<>();
            // get the all references of the reference Set/module
            MongoCursor<Document> genotypingDataCursor = MongoTemplateManager.get(id).getCollection(MongoTemplateManager.getMongoCollectionName(Sequence.class)).find().iterator();
            Document seq;

            String concatId = "";

            if (genotypingDataCursor != null) {

                List<String> listChecksum = new ArrayList<>();

                while (genotypingDataCursor.hasNext()) {
                    seq = genotypingDataCursor.next();
                    listChecksum.add((String) seq.get(Sequence.FIELDNAME_CHECKSUM));
                }
                // sort in lexicographic order
                Collections.sort(listChecksum);
                for (String checksum : listChecksum) {
                    concatId = concatId + checksum;
                }
                genotypingDataCursor.close();
            }

            String taxon = MongoTemplateManager.getTaxonName(id);
            String species = MongoTemplateManager.getSpecies(id);
            String taxoDesc = (species != null ? "Species: " + species : "") + (taxon != null && !taxon.equals(species) ? (species != null ? " ; " : "") + "Taxon: " + taxon : "");
            String refCountDesc;
            List<Assembly> assemblies = mongoTemplate.findAll(Assembly.class);
            List<String> assemblyNames = new ArrayList<>();
            MongoCollection<Document> projectColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(GenotypingProject.class));
//            if (assemblies.isEmpty())
//                    refCountDesc = projectColl.distinct(GenotypingProject.FIELDNAME_SEQUENCES, String.class).into(new ArrayList<>()).size() + " references ; ";
//            else {
                    refCountDesc = "";
                    for (Assembly assembly : assemblies) {
                        refCountDesc += (refCountDesc.isEmpty() ? "" : ", ") + projectColl.distinct(Assembly.getThreadBoundProjectContigsPath(), String.class).into(new ArrayList<>()).size() + " references (assembly " + assembly.getName() + ")";
                        assemblyNames.add(assembly.getName());
                    }
                    refCountDesc += " ; ";
//            }
            
            referenceSet = ReferenceSet.newBuilder()
                .setId(id)
                .setName(id)
                .setMd5checksum(Helper.convertToMD5(concatId))
                .setSourceAccessions(list)
                .setNcbiTaxonId(MongoTemplateManager.getTaxonId(id))
                .setAssemblyId(assemblyNames.isEmpty() ? null : StringUtils.join(assemblyNames, ", "))
                .setDescription((taxoDesc.isEmpty() ? "" : (taxoDesc + " ; ")) + refCountDesc + mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class)).estimatedDocumentCount() + " markers")
                .build();
        }
        return referenceSet;
    }

    @Override
    public Reference getReference(String id) throws AvroRemoteException {

        Reference reference = null;

        // get information from id
        String[] info = Helper.getInfoFromId(id, 3);
        if (info == null) {

            // wrong number of param or wrong module name
        } else {
            String module = info[0];
            int projId = Integer.parseInt(info[1]);
            String name = info[2];

            MongoTemplate mongoTemplate = MongoTemplateManager.get(module);

            GenotypingProject proj = mongoTemplate.findById(projId, GenotypingProject.class);
            Set<String> listRef = proj.getContigs(Assembly.getThreadBoundAssembly());

            // check if the sequence is in the list
            if (listRef.contains(name)) {

                long length = 0L;
                String checksum = Helper.convertToMD5("");

                Sequence sequence = mongoTemplate.findById(name, Sequence.class);
                if (sequence != null) {
                    length = sequence.getLength();
                    checksum = sequence.getChecksum();
                }

                List<String> liste = new ArrayList<>();
                // Checksum : MD5 of the upper-case sequence excluding all whitespace characters
                // length == 0 since we don't have this information in VCF files
                reference = Reference.newBuilder().setId(id)
                        .setMd5checksum(checksum)
                        .setName(name)
                        .setLength(length)
                        .setSourceAccessions(liste)
                        .build();
            }
        }
        return reference;
    }

    @Override
    public ListReferenceBasesResponse getReferenceBases(String id, ListReferenceBasesRequest lrbr) throws AvroRemoteException
    {
           String[] info = Helper.getInfoFromId(id, 3);
           String module = info[0];
           String seqName = info[2];
           String sequenceBases = "";
           int spaceCode = (int) '\n';
           int chevCode = (int) '>';

           if (lrbr.getEnd() > lrbr.getStart()) {
               Sequence seq = MongoTemplateManager.get(module).findById(seqName, Sequence.class);
               BufferedReader br = null;
               StringBuilder builder;
               if ((seq != null)) {
                   try {
                       br = new BufferedReader(new InputStreamReader(new FileInputStream(seq.getFilePath())));
                       builder = new StringBuilder();
                       String line;
                       int c;
                       int pos = 0;
                       // skip line to reach sequence
                   while ((line = br.readLine()) != null && !line.startsWith(">" + id)) {

                   }
                   // skip char to reach start pos
                   while ((c = br.read()) != -1 && pos < lrbr.getStart()) {
                       pos++;
                   }
                   builder.append((char) c);

                   while ((c = br.read()) != -1 && pos < lrbr.getEnd() && c != chevCode) {
                       if (c != spaceCode) {
                           builder.append((char) c);
                           pos++;
                       }
                   }
                   sequenceBases = builder.toString();

               }
               catch (IOException ex)
               {
                   LOG.warn("could not open file : " + ex);
               }
               finally
               {
                   try
                   {
                       if (br != null)
                           br.close();
                   }
                   catch (IOException ex)
                   {
                       LOG.warn("could not close writer : " + ex);
                   }
               }
           }
           else
               throw new GAException("No fasta for sequence " + seqName);
           }
           ListReferenceBasesResponse result = new ListReferenceBasesResponse();
           result.setSequence(sequenceBases);
           result.setOffset(lrbr.getStart());
           return result;
    }

    @Override
    public SearchCallSetsResponse searchCallSets(SearchCallSetsRequest scsr) throws AvroRemoteException {
        // get information from id
        String[] info = Helper.getInfoFromId(scsr.getVariantSetId(), 2);
        if (info == null)
            return null;

        GigwaSearchCallSetsRequest gscsr = (GigwaSearchCallSetsRequest) scsr;
        Authentication auth = tokenManager.getAuthenticationFromToken(tokenManager.readToken(gscsr.getRequest()));

        CallSet callSet;
        int start;
        int end;
        int pageSize;
        int pageToken = 0;
        String nextPageToken;

        String module = info[0];

        LinkedHashMap<String, Individual> indMap = mgdbDao.loadIndividualsWithAllMetadata(module, auth != null && auth.getAuthorities().contains(new SimpleGrantedAuthority(IRoleDefinition.ROLE_ADMIN)) ? null : AbstractTokenManager.getUserNameFromAuthentication(auth), Arrays.asList(Integer.parseInt(info[1])), null, null);

        List<CallSet> listCallSet = new ArrayList<>();
        int size = indMap.size();
        // if no pageSize specified, return all results
        if (scsr.getPageSize() != null) {
            pageSize = scsr.getPageSize();
        } else {
            pageSize = size;
        }
        if (scsr.getPageToken() != null) {
            pageToken = Integer.parseInt(scsr.getPageToken());
        }

        start = pageSize * pageToken;
        if (size - start <= pageSize) {
            end = size;
            nextPageToken = null;
        } else {
            end = pageSize * (pageToken + 1);
            nextPageToken = Integer.toString(pageToken + 1);
        }

        // create a callSet for each item in the list
        List<String> indList = new ArrayList() {{ addAll(indMap.keySet()); }};
        for (int i = start; i < end; i++) {
            final Individual ind = indMap.get(indList.get(i));
            CallSet.Builder csb = CallSet.newBuilder().setId(Helper.createId(module, info[1], ind.getId())).setName(ind.getId()).setVariantSetIds(Arrays.asList(scsr.getVariantSetId())).setSampleId(Helper.createId(module, info[1], ind.getId(), ind.getId()));

            if (!ind.getAdditionalInfo().isEmpty()) {
                            Map<String, String> addInfoMap = new HashMap();
                            for (String key:ind.getAdditionalInfo().keySet()) {
                                Object value = ind.getAdditionalInfo().get(key);
                                if (value instanceof String) {
                                    int spaces = ((String) value).length() - ((String) value).replaceAll(" ", "").length();
                                    if (spaces <= 5)
                                        addInfoMap.put(key, value.toString());
                                }
                            }
                            csb.setInfo(addInfoMap.keySet().stream().collect(Collectors.toMap(k -> k, k -> (List<String>) Arrays.asList(addInfoMap.get(k).toString()), (u,v) -> { throw new IllegalStateException(String.format("Duplicate key %s", u)); }, LinkedHashMap::new)));
                        }
            callSet = csb.build();
            listCallSet.add(callSet);
        }
        return SearchCallSetsResponse.newBuilder().setCallSets(listCallSet).setNextPageToken(nextPageToken).build();
    }

    @Override
    public SearchReferenceSetsResponse searchReferenceSets(SearchReferenceSetsRequest srsr) throws AvroRemoteException {
//    	long before = System.currentTimeMillis();
        List<String> list = new ArrayList<>();

        List<String> listModules = new ArrayList<>(MongoTemplateManager.getAvailableModules());
        Collections.sort(listModules);
        int start;
        int end;
        int pageSize;
        int pageToken = 0;
        String nextPageToken;

        int size = listModules.size();
        // if page size is not specified, return all results
        if (srsr.getPageSize() != null) {
            pageSize = srsr.getPageSize();
        } else {
            pageSize = size;
        }
        if (srsr.getPageToken() != null) {
            pageToken = Integer.parseInt(srsr.getPageToken());
        }

        start = pageSize * pageToken;

        // nextPageToken = null if no more result
        if (size - start <= pageSize) {
            end = size;
            nextPageToken = null;
        } else {
            end = pageSize * (pageToken + 1);
            nextPageToken = Integer.toString(pageToken + 1);
        }

        List<ReferenceSet> listRef = Collections.synchronizedList(new ArrayList<>());
        List<String> modulesToReturn = start != 0 || end != listModules.size() ? listModules.subList(start, end) : listModules;
	    if (!modulesToReturn.isEmpty()) {
	        List<ArrayList<String>> splitModuleCollections = Helper.evenlySplitCollection(modulesToReturn, Runtime.getRuntime().availableProcessors() * 2);
	        ExecutorService executor = Executors.newFixedThreadPool(splitModuleCollections.size());
	        for (int i=0; i<splitModuleCollections.size(); i++) {
	            Collection<String> modules = splitModuleCollections.get(i);
	            Thread t = new Thread() {
	                public void run() {
	                    // add a Reference Set for each existing module
	                    for (String module : modules) {
	                        MongoTemplate mongoTemplate = MongoTemplateManager.get(module);
	                        String taxon = MongoTemplateManager.getTaxonName(module);
	                        String species = MongoTemplateManager.getSpecies(module);
	                        String taxoDesc = (species != null ? "Species: " + species : "") + (taxon != null && !taxon.equals(species) ? (species != null ? " ; " : "") + "Taxon: " + taxon : "");
	                        String refCountDesc;
	                        List<Assembly> assemblies = mongoTemplate.findAll(Assembly.class);
	                        List<String> assemblyNames = new ArrayList<>();
	                        MongoCollection<Document> projectColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(GenotypingProject.class));
	                        if (assemblies.isEmpty())
	                                refCountDesc = projectColl.distinct(GenotypingProject.FIELDNAME_SEQUENCES, String.class).into(new ArrayList<>()).size() + " references ; ";
	                        else {
	                                refCountDesc = "";
	                                for (Assembly assembly : assemblies) {
	                                    refCountDesc += (refCountDesc.isEmpty() ? "" : ", ") + projectColl.distinct(GenotypingProject.FIELDNAME_CONTIGS + "." + assembly.getId(), String.class).into(new ArrayList<>()).size() + " references (assembly " + assembly.getName() + ")";
	                                    assemblyNames.add(assembly.getName());
	                                }
	                                refCountDesc += " ; ";
	                        }
	                        ReferenceSet referenceSet = ReferenceSet.newBuilder()
	                            .setId(module)
	                            .setName(module)
	                            .setMd5checksum("")    /* not supported at the time */
	                            .setSourceAccessions(list)
	                            .setNcbiTaxonId(MongoTemplateManager.getTaxonId(module))
	                            .setDescription((taxoDesc.isEmpty() ? "" : (taxoDesc + " ; ")) + refCountDesc + mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class)).estimatedDocumentCount() + " markers")
	                            .build();
	                        listRef.add(referenceSet);
	                    }
	                }
	            };
	            executor.execute(t);
	        }
	        executor.shutdown();
	        try {
	            executor.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
	        } catch (InterruptedException e) {
	            throw new AvroRemoteException(e);
	        }
	    }

//	    LOG.debug("searchReferenceSets took " + (System.currentTimeMillis() - before) + "ms");
        return SearchReferenceSetsResponse.newBuilder().setReferenceSets(listRef).setNextPageToken(nextPageToken).build();
    }

    @Override
    public SearchVariantSetsResponse searchVariantSets(SearchVariantSetsRequest svsr) throws AvroRemoteException {

        SearchVariantSetsResponse response = null;

        String[] info = Helper.getInfoFromId(svsr.getDatasetId(), 1);
        if (info == null)
            return null;

        String module = info[0];
        int start;
        int end;
        int pageSize;
        int pageToken = 0;

        String nextPageToken;

        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_NAME);
        q.fields().include(GenotypingProject.FIELDNAME_DESCRIPTION);
        q.fields().include(GenotypingProject.FIELDNAME_TECHNOLOGY);
        List<GenotypingProject> listProj = MongoTemplateManager.get(module).find(q, GenotypingProject.class);
        List<VariantSet> listVariantSet = new ArrayList<>();

        int size = listProj.size();
        // if page size is not specified, return all results
        if (svsr.getPageSize() != null) {
            pageSize = svsr.getPageSize();
        } else {
            pageSize = size;
        }
        if (svsr.getPageToken() != null) {
            pageToken = Integer.parseInt(svsr.getPageToken());
        }

        start = pageSize * pageToken;
        if (size - start <= pageSize) {
            end = size;
            nextPageToken = null;
        } else {
            end = pageSize * (pageToken + 1);
            nextPageToken = Integer.toString(pageToken + 1);
        }

        for (int i = start; i < end; i++) {
            GenotypingProject proj = listProj.get(i);
            String projId = Integer.toString(proj.getId());
            List<VariantSetMetadata> metadata = getMetadataList(module, projId);
            if (proj.getDescription() != null) {
                VariantSetMetadata vsmd = new VariantSetMetadata();
                vsmd.setKey(AbstractVariantData.VCF_CONSTANT_DESCRIPTION);
                vsmd.setValue(proj.getDescription());
                metadata.add(vsmd);
            }
            if (proj.getTechnology() != null && !proj.getTechnology().isEmpty()) {
                VariantSetMetadata vsmd = new VariantSetMetadata();
                vsmd.setKey(Constants.GENOTYPING_TECHNOLOGY);
                vsmd.setValue(proj.getTechnology());
                metadata.add(vsmd);
            }
            VariantSet variantSet = VariantSet.newBuilder()
                .setId(Helper.createId(module, projId))
                .setReferenceSetId(module)
                .setDatasetId(module)
                .setName(proj.getName())
                .setMetadata(metadata) // get the metadata from vcf header
                .build();
            listVariantSet.add(variantSet);
        }
        response = SearchVariantSetsResponse.newBuilder()
                .setVariantSets(listVariantSet)
                .setNextPageToken(nextPageToken)
                .build();

        return response;
    }

    @Override
    public GigwaSearchVariantsResponse searchVariants(SearchVariantsRequest svr) throws AvroRemoteException {

        GigwaSearchVariantsResponse response = null;
        // get extra info
        MgdbSearchVariantsRequest gsvr = (MgdbSearchVariantsRequest) svr;
        String info[] = Helper.getInfoFromId(svr.getVariantSetId(), 2);
        if (info == null) {
            // wrong number of param or wrong module name
        } else {
            boolean doCount = false;
            boolean doSearch = false;
            boolean doBrowse = false;
            boolean getGT = gsvr.isGetGT();
            // always do count because we need it for pagination mechanism.
            // As count are stored by query hash, it should not be a problem
            switch (gsvr.getSearchMode()) {
                case 0:
                    doCount = true;
                    doSearch = false;
                    doBrowse = false;
                    break;
                case 1:
                    doCount = false;
                    doSearch = true;
                    doBrowse = true;
                    break;
                case 2:
                    doCount = false;
                    doSearch = false;
                    doBrowse = true;
                    break;
                case 3:
                    doCount = true;
                    doSearch = true;
                    doBrowse = true;
                    break;
            }
            String module = info[0];
            int projId = Integer.parseInt(info[1]);

            Long count = null;
            long globalCount;
            
            String refPosPath = Assembly.getThreadBoundVariantRefPosPath(), pjContigsPath = Assembly.getThreadBoundProjectContigsPath();

            MongoCursor<Document> cursor = null;
            try
            {
                MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);
//                mongoTemplate.getCollection(mongoTemplate.getCollectionName(CachedCount.class)).drop();

                if (doSearch) {
                    // create a temp collection to store the result of the request
                    count = findVariants(gsvr);
                }
                else if (doCount || doBrowse) {
                    count = countVariants(gsvr, doBrowse);
                }

                if (count > 0 && doBrowse) {
                    String token = tokenManager.readToken(gsvr.getRequest());

                    if (token == null)
                        return null;

                    MongoCollection<Document> tempVarColl = MongoTemplateManager.getTemporaryVariantCollection(info[0], token, false, false, false);
                    FindIterable<Document> iterable;
                    if (gsvr.getSelectedVariantIds() != null && tempVarColl.countDocuments() > 0) {
                        iterable = tempVarColl.find(); //when searching on variant IDs, retrieving all temporary collection
                    } else {
                        VariantQueryWrapper varQueryWrapper = VariantQueryBuilder.buildVariantDataQuery(gsvr/*, getSequenceIDsBeingFilteredOn(gsvr.getRequest().getSession(), info[0])*/, true);
                        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();

                        //in this case, there is only one variantQueryDBList (no filtering on variant ids)
                        BasicDBList variantQueryDBList = !variantDataQueries.isEmpty() ? variantDataQueries.iterator().next() : new BasicDBList();

                        MongoCollection<Document> varCollForBuildingRows = tempVarColl.countDocuments() == 0 ? mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class)) : tempVarColl;
                        iterable = varCollForBuildingRows.find(!variantQueryDBList.isEmpty() ? new BasicDBObject("$and", variantQueryDBList) : new BasicDBObject());
                    }
                    iterable.collation(IExportHandler.collationObj);

                    if (gsvr.getSortBy() != null && gsvr.getSortBy().length() > 0)
                        iterable.sort(new BasicDBObject((!"_id".equals(gsvr.getSortBy()) ? refPosPath + "." : "") + gsvr.getSortBy(), Integer.valueOf("DESC".equalsIgnoreCase(gsvr.getSortDir()) ? -1 : 1)));
                    else if (mongoTemplate.findOne(new Query(new Criteria().andOperator(Criteria.where("_id").is(projId), Criteria.where(pjContigsPath + ".0").exists(true))), GenotypingProject.class) != null)
                        iterable.sort(new Document(refPosPath + "." + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPath + "." + ReferencePosition.FIELDNAME_START_SITE, 1));
                    else
                        iterable.sort(new Document("_id", 1));  // no positions available in this project: let's sort variants by ID

                    iterable.skip(Integer.parseInt(gsvr.getPageToken()) * gsvr.getPageSize()).limit(gsvr.getPageSize());    // skip the results we don't want

                    cursor = iterable.iterator();
                }
            }
            catch (Exception ex)
            {
                LOG.error("Error searching variants", ex);
                throw new GAException(ex);
            }
            globalCount = count == null ? 0 : count;
            // get the cursor containing variant from previously created temp collection
            // return null if cursor is empty
            if (cursor != null && cursor.hasNext()) {
                // we need to get callSet name and position in the callSet list to get corresponding genotype
                // if we don't want to retrieve genotype, just send an empty individuals list?
                Collection<GenotypingSample> samples = new ArrayList<>();
                if (getGT) {
                    try {
                        samples = MgdbDao.getSamplesForProject(module, projId, gsvr.getCallSetIds().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toList()));
                    } catch (ObjectNotFoundException ex) {
                        java.util.logging.Logger.getLogger(GigwaGa4ghServiceImpl.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                List<Variant> listVar = getVariantListFromDBCursor(module, Integer.parseInt(info[1]), cursor, samples);
                String nextPageToken = null;

                // if there is still more result after PageSize iterations
                if (globalCount > gsvr.getPageSize() * (Integer.parseInt(gsvr.getPageToken()) + 1)) {
                    nextPageToken = Integer.toString(Integer.parseInt(gsvr.getPageToken()) + 1);
                }
                response = new GigwaSearchVariantsResponse();
                response.setNextPageToken(nextPageToken);
                response.setVariants(listVar);
                    if (gsvr.getSearchMode() == 3) {
                        response.setCount(count);
                    }
                cursor.close();
                } else {
                    response = new GigwaSearchVariantsResponse();
                    response.setNextPageToken(null);
                    response.setVariants(new ArrayList<>());
                    response.setCount(count);
                }
            }
        return response;
    }

    @Override
    public SearchReferencesResponse searchReferences(SearchReferencesRequest srr) throws AvroRemoteException {

        SearchReferencesResponse response = null;

        // get information from id
        String[] info = Helper.getInfoFromId(srr.getReferenceSetId(), 1);
        if (info == null) {

            // wrong number of param or wrong module name
        } else {
            String module = info[0];
            GigwaSearchReferencesRequest gsr = (GigwaSearchReferencesRequest) srr;
            int projId = -1; // default : return all sequences of a module
            String[] array = gsr.getVariantSetId().split(Helper.ID_SEPARATOR);
            if (array.length > 1) {
                projId = Integer.parseInt(gsr.getVariantSetId().split(Helper.ID_SEPARATOR)[1]);
            }

            MongoTemplate mongoTemplate = MongoTemplateManager.get(module);
            int start;
            int end;
            int pageSize;
            int pageToken = 0;
            String nextPageToken;

            List<Reference> listReference = new ArrayList<>();
            List<String> accessions = new ArrayList<>();
            Map<String, Integer> mapSeq = new TreeMap<>(new AlphaNumericComparator());

            // allow search on checksum
            // but here only on one checksum v0.6.1 ?
            if (gsr.getMd5checksum() != null) {
                BasicDBObject query = new BasicDBObject();
                query.put(Sequence.FIELDNAME_CHECKSUM, gsr.getMd5checksum());
//                Document seq = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(Sequence.class)).find(query).first();
                // listSequence.add((String) seq.get("_id"));
            } else {
                Query q = new Query();
//                q.fields().include(GenotypingProject.FIELDNAME_SEQUENCES);
                if (projId != -1)
                    q.addCriteria(Criteria.where("_id").is(projId));
                List<GenotypingProject> listProj = mongoTemplate.find(q, GenotypingProject.class);
                for (int i = 0; i < listProj.size(); i++) {
                    for (String seq : listProj.get(i).getContigs(Assembly.getThreadBoundAssembly())) {
                        mapSeq.put(seq, i + 1);
                    }
                }
            }
            int size = mapSeq.size();
            // if page size is not specified, return all results
            if (gsr.getPageSize() != null) {
                pageSize = gsr.getPageSize();
            } else {
                pageSize = size;
            }
            if (gsr.getPageToken() != null) {
                pageToken = Integer.parseInt(gsr.getPageToken());
            }
            start = pageSize * pageToken;
            if (size - start <= pageSize) {
                end = size;
                nextPageToken = null;
            } else {
                end = pageSize * (pageToken + 1);
                nextPageToken = Integer.toString(pageToken + 1);
            }
            ArrayList<BasicDBObject> pipeline = new ArrayList<>();
            pipeline.add(new BasicDBObject("$match", new BasicDBObject("_id", new BasicDBObject("$in", mapSeq.keySet()))));
            MongoCursor<Document> sqCursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(Sequence.class)).aggregate(pipeline).allowDiskUse(true).iterator();

            Iterator<String> iteratorName = mapSeq.keySet().iterator();
            Iterator<Integer> iteratorId = mapSeq.values().iterator();

            // create and add the corresponding Reference for each sequence
            for (int i = start; i < end; i++) {
                Document sequence = !sqCursor.hasNext() ? null : sqCursor.next();
                String name = iteratorName.next();
                String projectId = Integer.toString(iteratorId.next());

                // Checksum : MD5 of the upper-case sequence excluding all whitespace characters (we usually don't have it)
                String md5 = sequence == null ? Helper.convertToMD5("") : (String) sequence.get(Sequence.FIELDNAME_CHECKSUM);
                Reference reference = Reference.newBuilder().setId(Helper.createId(module, projectId, name))
                        .setMd5checksum(md5 == null ? Helper.convertToMD5("") : md5)
                        .setName(name)
                        .setLength(sequence == null ? 0 : (long) sequence.get(Sequence.FIELDNAME_LENGTH))    // length == 0 when we don't have this information
                        .setSourceAccessions(accessions)
                        .build();

                listReference.add(reference);
            }

            response = SearchReferencesResponse.newBuilder()
                    .setReferences(listReference)
                    .setNextPageToken(nextPageToken)
                    .build();
        }
        return response;
    }

    /**
     * return the ID of an ontology term
     *
     * @param name
     * @return
     */
    public String getOntologyId(String name) {

        return MongoTemplateManager.getOntologyMap().get(name);
    }

    /**
     * Get annotations for a specific variant Waiting for ga4gh schema work only
     * for VCF 4.2 version using ANN notation
     *
     * @param id variant ID
     * @return snpEff annotation for this variant
     */
    public VariantAnnotation getVariantAnnotation(String id) {
        VariantAnnotation.Builder variantAnnotationBuilder = VariantAnnotation.newBuilder()
            .setVariantId(id)
            .setId(id)
            .setVariantAnnotationSetId(id.substring(0, id.lastIndexOf(Helper.ID_SEPARATOR))); // which variant annotation set?

        // get information from id
        String[] info = Helper.getInfoFromId(id, 3);
        if (info == null) {
            // wrong number of param or wrong module name
        } else {
            String module = info[0];
            String variantId = info[2];
            String[] headerField;
            String header;

            // parse annotation fields
            BasicDBObject queryVarAnn = new BasicDBObject();
            BasicDBObject varAnnField = new BasicDBObject();
            queryVarAnn.put("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, variantId);
            varAnnField.put(VariantData.FIELDNAME_KNOWN_ALLELES, 1);
            varAnnField.put(VariantData.SECTION_ADDITIONAL_INFO, 1);
            Document variantRunDataObj = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).find(queryVarAnn).projection(varAnnField)
                .sort(new BasicDBObject(AbstractVariantData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME, -1))  /*FIXME: this method should be called separately for each run*/
                .first();
            Document variantAnnotationObj = variantRunDataObj != null ? (Document) variantRunDataObj.get(VariantRunData.SECTION_ADDITIONAL_INFO) : null;
            if (variantAnnotationObj != null)
            {
                String ann = (String) variantAnnotationObj.get(VcfImport.ANNOTATION_FIELDNAME_ANN);
                if (ann == null)
                    ann = (String) variantAnnotationObj.get(VcfImport.ANNOTATION_FIELDNAME_CSQ);
                boolean fAnnStyle = ann != null;
                if (!fAnnStyle)
                    ann = (String) variantAnnotationObj.get(VcfImport.ANNOTATION_FIELDNAME_EFF);
                Map<String, List<String>> additionalInfo = new HashMap<>();

                String[] tableTranscriptEffect = new String[0];
                if (ann != null)
                {
                    tableTranscriptEffect = ann.split(",");
                    List<TranscriptEffect> transcriptEffectList = new ArrayList<>();

                    // source version is stored in the ontology map
                    String sourceVersion = getOntologyId(Constants.VERSION) == null ? "" : getOntologyId(Constants.VERSION);

                    BasicDBObject fieldHeader = new BasicDBObject(AbstractVariantData.VCF_CONSTANT_INFO_META_DATA + "." + (fAnnStyle ? VcfImport.ANNOTATION_FIELDNAME_ANN : VcfImport.ANNOTATION_FIELDNAME_EFF) + "." + AbstractVariantData.VCF_CONSTANT_DESCRIPTION, 1);
                    if (fAnnStyle)
                        fieldHeader.put(AbstractVariantData.VCF_CONSTANT_INFO_META_DATA + "." + VcfImport.ANNOTATION_FIELDNAME_CSQ + "." + AbstractVariantData.VCF_CONSTANT_DESCRIPTION, 1);

                    MongoCollection<Document> vcfHeaderColl = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class));
                    BasicDBList vcfHeaderQueryOrList = new BasicDBList();
                    for (String key : fieldHeader.keySet())
                        vcfHeaderQueryOrList.add(new BasicDBObject(key, new BasicDBObject("$exists", true)));

                    Document vcfHeaderEff = vcfHeaderColl.find(new BasicDBObject("$or", vcfHeaderQueryOrList)).projection(fieldHeader).first();
                    if (vcfHeaderEff != null) {
	                    ArrayList<String> headerList = new ArrayList<>();
	                    LinkedHashSet<String> usedHeaderSet = new LinkedHashSet<>();
	                    if (!fAnnStyle)
	                        headerList.add("Consequence");    // EFF style annotations
	                    Document annInfo = (Document) ((Document) vcfHeaderEff.get(AbstractVariantData.VCF_CONSTANT_INFO_META_DATA)).get(fAnnStyle ? VcfImport.ANNOTATION_FIELDNAME_ANN : VcfImport.ANNOTATION_FIELDNAME_EFF);
	                    if (annInfo == null && fAnnStyle)
	                        annInfo = (Document) ((Document) vcfHeaderEff.get(AbstractVariantData.VCF_CONSTANT_INFO_META_DATA)).get(VcfImport.ANNOTATION_FIELDNAME_CSQ);
	                    if (annInfo != null) {
	                        header = (String) annInfo.get(AbstractVariantData.VCF_CONSTANT_DESCRIPTION);
	                        if (header != null) {
	                            // consider using the headers for additional info keySet
	                            String sBeforeFieldList = fAnnStyle ? ": " : " (";
	                            headerField = header.substring(header.indexOf(sBeforeFieldList) + sBeforeFieldList.length(), fAnnStyle ? header.length() : header.indexOf(")")).replaceAll("'", "").split("\\|");
	                            for (String head : headerField)
	                                headerList.add(head.replace("[", "").replace("]", "").trim());
	                        }
	                    }
	
	                    List<AnalysisResult> listAnalysisResults = new ArrayList<>();
	
	                    for (int i=0; i<tableTranscriptEffect.length; i++) {
	                        ArrayList<String> values = new ArrayList<>();
	
	                        if (!fAnnStyle) {    // EFF style annotations
	                            int parenthesisPos = tableTranscriptEffect[i].indexOf("(");
	                            values.add(tableTranscriptEffect[i].substring(0, parenthesisPos));
	                            tableTranscriptEffect[i] = tableTranscriptEffect[i].substring(parenthesisPos + 1).replace(")", "");
	                        }
	
	                        List<OntologyTerm> ontologyList = new ArrayList<>();
	
	                        String[] effectFields = tableTranscriptEffect[i].split("\\|", -1);
	                        for (int j=0; j<effectFields.length; j++)
	                        {
	                            values.add(effectFields[j]);
	                            if (effectFields[j].endsWith(")"))
	                            {
	                                String[] splitVal = effectFields[j].substring(0,  effectFields[j].length() - 1).split("\\(");
	                                if (splitVal.length == 2)
	                                    try
	                                    {
	                                        AnalysisResult analysisResult = new AnalysisResult();
	                                        analysisResult.setAnalysisId(headerList.get(j));
	                                        analysisResult.setResult(splitVal[0]);
	                                        analysisResult.setScore((int)(100 * Float.parseFloat(splitVal[1])));
	                                        listAnalysisResults.add(analysisResult);
	                                    }
	                                    catch (NumberFormatException ignored)
	                                    {}
	                            }
	                        }
	
	                        int impactIndex = headerList.indexOf(fAnnStyle ? "IMPACT" : "Effefct_Impact");
	                        if (impactIndex != -1)
	                        {
	                            String[] impact = values.get(impactIndex).split("&");
	
	                            for (String anImpact : impact) {
	                                String ontologyId = getOntologyId(anImpact);
	                                if (ontologyId == null) {
	                                    ontologyId = "";
	                                }
	                                OntologyTerm ontologyTerm = OntologyTerm.newBuilder()
	                                        .setId(ontologyId)
	                                        .setSourceName("sequence ontology")
	                                        .setSourceVersion(sourceVersion)
	                                        .setTerm(anImpact)
	                                        .build();
	                                ontologyList.add(ontologyTerm);
	                            }
	                        }
	
	                        HGVSAnnotation.Builder hgvsBuilder = HGVSAnnotation.newBuilder();
	                        AlleleLocation cDnaLocation = null;
	                        AlleleLocation cdsLocation = null;
	                        AlleleLocation proteinLocation = null;
	                        int nC = headerList.indexOf("HGVSc"), nP = headerList.indexOf("HGVSp"), nT = headerList.indexOf("Transcript");
	                        if ((nC != -1 && !values.get(nC).isEmpty()) || (nP != -1 && !values.get(nP).isEmpty()) || (nT != -1 && !values.get(nT).isEmpty()))
	                        {
	                            if (nC != -1)
	                                hgvsBuilder.setGenomic(values.get(nC));
	                            if (nT != -1)
	                                hgvsBuilder.setTranscript(values.get(nT));
	                            if (nP != -1)
	                                hgvsBuilder.setProtein(values.get(nP));
	                        }
	
	                        if (fAnnStyle)
	                        {
	                            for (String positionHeader : Arrays.asList("cDNA_position", "CDS_position", "Protein_position"))
	                            {
	                                int nPos = headerList.indexOf(positionHeader);
	                                if (nPos != -1)
	                                {
	                                    String value = values.get(nPos);
	                                    if (!value.equals(""))
	                                    {
	                                        AlleleLocation.Builder allLocBuilder = AlleleLocation.newBuilder();
	                                        String[] splitVals = value.split("/");
	                                        if (splitVals.length == 1 && value.contains("-"))
	                                            splitVals = value.split("-");    // sometimes used as separator
	                                        try
	                                        {
	                                            allLocBuilder.setStart(Integer.parseInt(splitVals[0]));
	                                        }
	                                        catch (NumberFormatException ignored)
	                                        {}
	
	                                        if (allLocBuilder.getStart() == 0)
	                                            continue;
	
	                                        boolean fWorkingOnProtein = "Protein_position".equals(positionHeader);
	
	                                        String sRefAllele = ((List<String>) variantRunDataObj.get(VariantData.FIELDNAME_KNOWN_ALLELES)).get(0);
	                                        if (!fWorkingOnProtein)
	                                            allLocBuilder.setEnd(allLocBuilder.getStart() + sRefAllele.length() - 1);
	//                                        else
	                                            /* TODO: don't know how to calculate END field for proteins */
	
	                                        if ("cDNA_position".equals(positionHeader))
	                                            cDnaLocation = allLocBuilder.build();
	                                        else if (!fWorkingOnProtein)
	                                            cdsLocation = allLocBuilder.build();
	                                        else
	                                            proteinLocation = allLocBuilder.build();
	                                    }
	                                }
	                            }
	                        }
	
	                        TranscriptEffect transcriptEffect = TranscriptEffect.newBuilder()
	                                .setAlternateBases(values.get(0))
	                                .setId(id + Helper.ID_SEPARATOR + i)
	                                .setEffects(ontologyList)
	                                .setHgvsAnnotation(hgvsBuilder.build())
	                                .setCDNALocation(cDnaLocation)
	                                .setCDSLocation(cdsLocation)
	                                .setProteinLocation(proteinLocation)
	                                .setFeatureId(values.get(6))
	                                .setAnalysisResults(listAnalysisResults)
	                                .build();

	                        transcriptEffectList.add(transcriptEffect);
	                        additionalInfo.put(Constants.ANN_VALUE_LIST_PREFIX + i, values);
	                        for (int j=0; j<values.size(); j++)
	                            if (!values.get(j).isEmpty())
	                                usedHeaderSet.add(headerList.get(j));
	                    }
	
	                    for (int i=0; i<tableTranscriptEffect.length; i++)
	                    {
	                        List<String> keptValues = new ArrayList<String>(), allValues = additionalInfo.get(Constants.ANN_VALUE_LIST_PREFIX + i);
	                        for (int j=0; j<allValues.size(); j++)
	                            if (usedHeaderSet.contains(headerList.get(j)))
	                                keptValues.add(allValues.get(j));
	                        additionalInfo.put(Constants.ANN_VALUE_LIST_PREFIX + i, keptValues);
	                    }

	                    ArrayList<String> properlySortedUsedHeaderList = new ArrayList<>();
	                    for (String aHeader : headerList)
	                        if (usedHeaderSet.contains(aHeader))
	                            properlySortedUsedHeaderList.add(aHeader);
	                    additionalInfo.put(Constants.ANN_HEADER, new ArrayList<String>(properlySortedUsedHeaderList));

	                    variantAnnotationBuilder.setTranscriptEffects(transcriptEffectList);                    }
                }

                TreeMap<String, String> metadata = new TreeMap<>();
                for (String key : variantAnnotationObj.keySet())
                    // do not store EFF_ge / EFF_nm / EFF / ANN / CSW
                    if (!key.equals(VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_GENE) && !key.equals(VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME) && !key.equals(VcfImport.ANNOTATION_FIELDNAME_ANN) && !key.equals(VcfImport.ANNOTATION_FIELDNAME_CSQ) && !key.equals(VcfImport.ANNOTATION_FIELDNAME_EFF) && !key.equals(""))
                        metadata.put(key, variantAnnotationObj.get(key).toString());
                additionalInfo.put(Constants.METADATA_HEADER, new ArrayList<String>(metadata.keySet()));
                additionalInfo.put(Constants.METADATA_VALUE_LIST, new ArrayList<String>(metadata.values()));
                variantAnnotationBuilder.setInfo(additionalInfo);
            }
        }
        return variantAnnotationBuilder.build();
    }

    @Override
    public Map<String, String> getAnnotationHeaders(String module, int projId) {

        Map<String, String> annHeaders = new HashMap<>();
        BasicDBObject queryVarAnn = new BasicDBObject();
        BasicDBObject varAnnField = new BasicDBObject();
        queryVarAnn.put("_id." + DBVCFHeader.VcfHeaderId.FIELDNAME_PROJECT, projId);
        varAnnField.put(AbstractVariantData.VCF_CONSTANT_INFO_META_DATA, 1);
        varAnnField.put(AbstractVariantData.VCF_CONSTANT_INFO_FORMAT_META_DATA, 1);
        Document result = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class)).find(queryVarAnn).projection(varAnnField).first();
        if (result != null) {
            Document metaDataHeader = (Document) result.get(AbstractVariantData.VCF_CONSTANT_INFO_META_DATA);
            for (String key : metaDataHeader.keySet()) {
                annHeaders.put(key, (String) ((Document) metaDataHeader.get(key)).get(AbstractVariantData.VCF_CONSTANT_DESCRIPTION));
            }
            Document formatHeader = (Document) result.get(AbstractVariantData.VCF_CONSTANT_INFO_FORMAT_META_DATA);
            for (String key : formatHeader.keySet()) {
                annHeaders.put(key, (String) ((Document) formatHeader.get(key)).get(AbstractVariantData.VCF_CONSTANT_DESCRIPTION));
            }
        }
        return annHeaders;
    }

    @Override
    public TreeMap<String, HashMap<String, String>> getExportFormats() {
    	boolean enableExperimentalFeatures = Boolean.TRUE.equals(Boolean.parseBoolean(appConfig.get("enableExperimentalFeatures")));
        TreeMap<String, HashMap<String, String>> exportFormats = new TreeMap<>();
        try {
            for (IExportHandler exportHandler : Stream.of(AbstractIndividualOrientedExportHandler.getIndividualOrientedExportHandlers().values(), AbstractMarkerOrientedExportHandler.getMarkerOrientedExportHandlers().values()).flatMap(Collection::stream).collect(Collectors.toList()))
                if (enableExperimentalFeatures || !ExperimentalFeature.class.isAssignableFrom(exportHandler.getClass())) {
	            	HashMap<String, String> info = new HashMap<>();
	                info.put("desc", exportHandler.getExportFormatDescription());
	                info.put("supportedPloidyLevels", StringUtils.join(ArrayUtils.toObject(exportHandler.getSupportedPloidyLevels()), ";"));
	                info.put("dataFileExtensions", StringUtils.join(exportHandler.getExportDataFileExtensions(), ";"));
	                info.put("supportedVariantTypes", StringUtils.join(exportHandler.getSupportedVariantTypes(), ";"));
	                exportHandler.setTmpFolder(servletContext.getRealPath(FRONTEND_URL + File.separator + TMP_OUTPUT_EXTRACTION_FOLDER));
	                exportFormats.put(exportHandler.getExportFormatName(), info);
	            }
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException | SecurityException ex) {
            LOG.debug("error", ex);
        }
        return exportFormats;
    }

    public List<Comparable> searchVariantsLookup(String module, int projectId, String lookupText) throws AvroRemoteException {

        MongoTemplate mongoTemplate = MongoTemplateManager.get(module);

        List<Comparable> values = new ArrayList<>();

        MongoCollection<Document> collection = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantData.class));

        BasicDBObject whereQuery = new BasicDBObject();
        whereQuery.put("_id", Pattern.compile(".*\\Q" + lookupText + "\\E.*", Pattern.CASE_INSENSITIVE));
        whereQuery.put(VariantData.FIELDNAME_RUNS + "." + Run.FIELDNAME_PROJECT_ID, projectId);

        int maxSize = 50;
        try {
            String variantIdLookupMaxSize = appConfig.get("variantIdLookupMaxSize");
            maxSize = Integer.parseInt(variantIdLookupMaxSize);
        } catch (Exception e) {
            LOG.debug("can't read variantIdLookupMaxSize in config, using maxSize=50");
        }

        MongoCursor<Document> cursor = collection.aggregate(
            Arrays.asList(
                Aggregates.match(whereQuery),
                Aggregates.group("$_id"),
                Aggregates.limit(maxSize+1)
            )
        ).iterator();

        try {
            while (cursor.hasNext()) {
                values.add((Comparable) cursor.next().get("_id"));
            }
        } finally {
           cursor.close();
        }

        if (values.size() > maxSize)
            return Arrays.asList("Too many results, please refine search!");

        return values;
    }
    
    public List<String> searchGenesLookup(String module, int projectId, String lookupText) throws AvroRemoteException {
    	long before = System.currentTimeMillis();
    	String fieldPath = "_id";

        MongoTemplate mongoTemplate = MongoTemplateManager.get(module);
        List<String> values = new ArrayList<>();
        MongoCollection<Document> collection = mongoTemplate.getCollection(MgdbDao.COLLECTION_NAME_GENE_CACHE);

        BasicDBObject whereQuery = new BasicDBObject();
        whereQuery.put(fieldPath, Pattern.compile(".*" + lookupText + ".*", Pattern.CASE_INSENSITIVE));
        
        int maxSize = 50;
        try {
            String variantIdLookupMaxSize = appConfig.get("variantIdLookupMaxSize");
            maxSize = Integer.parseInt(variantIdLookupMaxSize);
        } catch (Exception e) {
            LOG.debug("Can't read variantIdLookupMaxSize in config, using maxSize=50");
        }
        
        MongoCursor<Document> cursor = collection.find(whereQuery).iterator();

            try {
                while (cursor.hasNext()) {
                    values.add((String) cursor.next().get("_id"));
                }
            } finally {
               cursor.close();
            }

            if (values.size() > maxSize)
                return Arrays.asList("Too many results, please refine search!");

//        DistinctIterable<String> distinctValues = collection.distinct(fieldPath, whereQuery, String.class);
//        
//        for (String value : distinctValues) {
//            values.add(value);
//            if (values.size() >= maxSize) {
//                break;
//            }
//        }
//
//        if (values.size() > maxSize) {
//            values.clear();
//            values.add("Too many results, please refine search!");
//        }

        LOG.info("searchGenesLookup found " + values.size() + " results in " + (System.currentTimeMillis() - before) / 1000d + "s");
        return values;
    }

	@Override
	public void setServletContext(ServletContext servletContext) {
		this.servletContext = servletContext;
	}
}
