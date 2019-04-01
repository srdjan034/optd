
void putPointerAfterColon(char ** s_tmp)
{
    while(**s_tmp != ':')
        (*s_tmp)++;
    (*s_tmp)++;
}

/*
    Lightweight Conf.json parser
 */
int readJsonIntValue(char * s, char * jsonFieldName, int * value)
{
    char * s_tmp = strstr(s, jsonFieldName);
    
    if(s_tmp != NULL)
    {
        int i;
        for(i = 0; i < strlen(jsonFieldName); i++)
            s_tmp++;
        
        putPointerAfterColon(&s_tmp);

        sscanf(s_tmp, "%d", value);
        
        return 1;
    }
    else
        return -1;
}

/*
    Lightweight Conf.json parser
 */
Conf * readConf(char fileName[])
{
    FILE * f = fopen(fileName, "rb");
    
    Conf * conf = (Conf *) malloc(sizeof(Conf)); 
    
    char * buffer = NULL;
    long length;
    
    if (f)
    {
        fseek (f, 0, SEEK_END);
        length = ftell (f);
        fseek (f, 0, SEEK_SET);
        buffer = malloc (length);
        
        if (buffer)
        {
          fread (buffer, 1, length, f);
        }
    }
    
    readJsonIntValue(buffer, "nodeCount", &conf->nodeCount);
    readJsonIntValue(buffer, "processorsPerNode", &conf->processorsPerNode);
    readJsonIntValue(buffer, "producerCount", &conf->producerCount);
    readJsonIntValue(buffer, "consumerCount", &conf->consumerCount);
    readJsonIntValue(buffer, "partialTrajectoryLength", &conf->partialTrajectoryLength);
    readJsonIntValue(buffer, "chunkSize", &conf->chunkSize);
    readJsonIntValue(buffer, "particleCount", &conf->particleCount);
    
    fclose(f);
    
    return conf;
}