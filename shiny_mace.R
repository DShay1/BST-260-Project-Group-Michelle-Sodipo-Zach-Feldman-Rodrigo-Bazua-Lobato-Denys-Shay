#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(ggthemes)
library(rstudioapi)
library(broom)
library(readr)
library(scales)
library(tidyverse)
library(tidyr)
library(splines2)
library(gam)
library(foreign)
library(Hmisc)
library(haven)
library(gtools)
library(splitstackshape)
library(ResourceSelection)
library(LogisticDx)
library(ROCR)
library(knitr)
library(nnet)
library(VGAM)
library(survival)
library(survminer)
library(gridExtra)
library(dplyr)
library(plyr)
library(splitstackshape)
library(caret)
library(e1071)
library(pROC)
library(MASS)
library(rpart)
library(randomForest)
library(betareg)
library(Hmisc)
library(haven)
library(stringr)
library(forcats)
library(viridis)
library(hrbrthemes)
library(tibble)
library(gdata)
library(reshape2)
library(shinythemes)
library(magrittr)
library(kableExtra)

load("dermptscomplete.RData")
load("admsubjcomplete.RData")

derm <- dermptscomplete
adm <- admsubjcomplete

data <- dermptscomplete %>% dplyr::select(subject_id,mace, gender, age_cat, ethnicity, mi, chf, hiv, pvd, cva, copd, rheum,  livmodsev, malig, miliv, pleg, dm, numadm) %>% 
    mutate(gender = as.factor(gender), 
           age_cat = as.factor(age_cat), 
           ethnicity = as.factor(ethnicity), 
           mi = as.factor(mi), 
           chf = as.factor(chf),
           hiv = as.factor(hiv), 
           pvd = as.factor(pvd), 
           cva = as.factor(cva), 
           copd = as.factor(copd), 
           rheum = as.factor(rheum), 
           livmodsev = as.factor(livmodsev),
           malig = as.factor(malig),
           miliv = as.factor(miliv),
           pleg = as.factor(pleg), 
           dm = as.factor(dm),
           numadm = as.numeric(numadm),
           mace = as.factor(mace)) 

# Use complete cases 
data <- data[complete.cases(data), ]

set.seed(1)

x <- stratified(data, "mace", 0.7, keep.rownames = TRUE)
train_set <- x %>% dplyr::select(-rn)
train_index <- as.numeric(x$rn)
test_set <- data[-train_index,]

###### KNN with k=4
knn_fit <- knn3(mace~.,data = dplyr::select(train_set, mace, gender, age_cat, ethnicity, mi, chf, hiv, pvd, cva, copd, rheum, livmodsev, malig, miliv, pleg, dm, numadm), k = 4)
f_hat_4 <- predict(knn_fit, newdata = test_set)[,2]
tab_4 <- table(pred=round(f_hat_4), truth=test_set$mace)
confusionMatrix(tab_4, positive = "1")

###### KNN with k=8
knn_fit <- knn3(mace~.,data = dplyr::select(train_set, mace, gender, age_cat, ethnicity, mi, chf, hiv, pvd, cva, copd, rheum, livmodsev, malig, miliv, pleg, dm, numadm), k = 8)
f_hat_8 <- predict(knn_fit, newdata = test_set)[,2]
tab_8 <- table(pred=round(f_hat_8), truth=test_set$mace)
confusion8 <- confusionMatrix(tab_8, positive = "1")

###### KNN with k=12
knn_fit <- knn3(mace~.,data = dplyr::select(train_set, mace, gender, age_cat, ethnicity, mi, chf, hiv, pvd, cva, copd, rheum, livmodsev, malig, miliv, pleg, dm, numadm), k =12)
f_hat_12 <- predict(knn_fit, newdata = test_set)[,2]
tab_12 <- table(pred=round(f_hat_12), truth=test_set$mace)
confusion12 <-confusionMatrix(tab_12, positive = "1")

#### LDA model
lda_fit <- lda(mace ~ gender + age_cat + ethnicity + mi + chf + hiv + pvd + cva + copd + rheum + livmodsev + malig + miliv + pleg + dm + numadm, data = train_set)
p_hat_lda <- predict(lda_fit, newdata = test_set)$posterior[,2]
y_hat_lda <- factor(ifelse(p_hat_lda > 0.5, 1, 0))
confusionMatrix(data = as.factor(y_hat_lda), reference = test_set$mace, positive = "1")

#### Random Forrest
RNGkind(sample.kind = "Rounding")
set.seed(1)
forest_fit <- randomForest(mace ~ gender + age_cat + ethnicity + mi + chf + hiv + pvd + cva + copd + rheum + livmodsev + malig + miliv + pleg + dm + numadm, data = train_set)
p_hat_forest <- predict(forest_fit, newdata = test_set, type = "prob")[, 2]
y_hat_forest <- factor(ifelse(p_hat_forest > 0.5, 1, 0))
confusionMatrix(data = as.factor(y_hat_forest), reference = test_set$mace, positive = "1")

#### ROC Curves
knn_4_roc <- roc(test_set$mace, f_hat_4)
knn_8_roc <- roc(test_set$mace, f_hat_8)
knn_12_roc <- roc(test_set$mace, f_hat_12)
roc_lda <- roc(test_set$mace, p_hat_lda)
roc_forest <- roc(test_set$mace, p_hat_forest)

####AUCs
aucknn4 <- auc(knn_4_roc)
aucknn8 <-     auc(knn_8_roc)
aucknn12 <-     auc(knn_12_roc)
aucknnlda <-     auc(roc_lda)
aucknnroc <-     auc(roc_forest)

glmmaceall <- glm(mace ~ age_cat + gender + ethnicity + language + admission_location + insurance + marital_status + hid + vit + psor + atop + alop + mi + chf + hiv + pvd + cva + copd + rheum + malig + miliv + livmodsev + dm, data = admsubjcomplete, family = "binomial")
tab <- cbind(round(exp(summary(glmmaceall)$coefficients[,1]),3),round(exp(confint.default(glmmaceall))[,1],3),round(exp(confint.default(glmmaceall))[,2],3),round((summary(glmmaceall)$coefficients[,4]),5))
colnames(tab) <- c("Odds ratio for MACE","2.5% CI limit","97.5% CI limit","p-value")
kable(tab)

glmmaceall <- glm(mace ~ age_cat + gender + ethnicity + language + admission_location + insurance + marital_status + hid + vit + psor + atop + alop + mi + chf + hiv + pvd + cva + copd + rheum + malig + miliv + livmodsev + dm, data = admsubjcomplete, family = "binomial")
tab2 <- cbind(round(exp(summary(glmmaceall)$coefficients[,1]),3),round(exp(confint.default(glmmaceall))[,1],3),round(exp(confint.default(glmmaceall))[,2],3),round((summary(glmmaceall)$coefficients[,4]),5))
colnames(tab2) <- c("Odds ratio for MACE","2.5% CI limit","97.5% CI limit","p-value")
kable(tab2)

glmmacederm <- glm(mace ~ age_cat + gender+ ethnicity + language + admission_location + insurance + marital_status + mi + chf + hiv + pvd + cva + copd + rheum + malig + miliv + livmodsev + dm, data = dermptscomplete, family = "binomial")
tab3 <- cbind(round(exp(summary(glmmacederm)$coefficients[,1]),3),round(exp(confint.default(glmmacederm))[,1],3),round(exp(confint.default(glmmacederm))[,2],3),round((summary(glmmacederm)$coefficients[,4]),5))
colnames(tab3) <- c("Odds ratio for MACE","2.5% CI limit","97.5% CI limit","p-value")
kable(tab3)


AUCS <- list("AUC Knn = 4" = aucknn4, "AUC Knn = 8" = aucknn8, "AUC Knn = 12" = aucknn12, "AUC LDA" = aucknnlda, "AUC RF" = aucknnroc )

# Define UI for application that draws a histogram
ui <- fluidPage(
    navbarPage("R Shiny Dashboard",
               tabPanel("Welcome",
                        tabName = "welcome",
                        icon=icon("door-open"),
                        
                        fluidPage(theme=shinytheme("cerulean"),
                                  h1("Welcome to our MACE prediction shiny dashboard"),
                                  br(),
                                  p(strong(tags$u("You can analyze the plots exploring the data"))),
                                  p("Choose between one of the plots"),  
                                  br(),
                                  p(strong(tags$u("You can review the regression models"))),
                                  p("Choose the regression model you want to review"),
                                  br(),
                                  p(strong(tags$u("You can compare machine learning models"))),
                                  p("Choose the model and review its Confusion Matrix")
                        )),
               tabPanel("EDA",
                        tabname="Tables and statistics",
                        icon=icon("bar-chart"),
                        headerPanel("Variables in the Study Population"), 

                        tabsetPanel(
                            tabPanel(
                                titlePanel(h4("Ethnicity",style='padding-left: 15px')),
                                tabPanel(title="Plots",
                                         h3("Race/Ethnicity of Study Population"),    
                                         plotOutput(outputId = "ethnicityplot"),
                                         br()               
                                )
                            )
                        
                            ,
                            tabPanel(
                                titlePanel(h4("Age by Gender", style='padding-left: 15px')),
                                tabPanel(title="Plots",
                                         h3("Age by Gender"),   
                                         plotOutput(outputId = "agegenderplot"),
                                         br()                
                                )
                            ),
                            tabPanel(
                                titlePanel(h4("Admission by Insurance",style='padding-left: 15px')),
                                tabPanel(title="Admission Type & Insurance",
                                         h3("Admission Type & Insurance"),    
                                         plotOutput(outputId = "adminsplot"),
                                         br()                
                                )
                            ),
                            tabPanel(
                                titlePanel(h4("MACE by ethnicity", style='padding-left: 15px')),
                                tabPanel(title="Plots",
                                         h3("MACE Outcomes by Age and Ethnicity"),    
                                         plotOutput(outputId = "maceethnicityplot"),
                                         br()                
                                )
                            )
                        )
                        ),
               tabPanel("Data Analysis",
                        tabname="regression",
                        icon=icon("calculator"),
                        headerPanel("Data Analysis"),
                        tabsetPanel(
                            tabPanel(
                                titlePanel(h4("CSD as a whole",style='padding-left: 15px')),
                                tabPanel(title="Plos",
                                         h3("Treating CISD conditions as a whole"),    
                                         tableOutput("dermwhole"),
                                         br(),
                                         p("Here, we can see that in addition to classical risk factors for MACE (such as age, male sex, prior MI, diagnosed CHF, and COPD), **having a CISD confers a 1.20-times higher odds of having a MACE during hospitalization**, holding the other covariates constant."),
                                         br()
                                )
                            ),
                            tabPanel(
                                titlePanel(h4("CISD individually",style='padding-left: 15px')),
                                tabPanel(title="Plos",
                                         h3("Treating CISD conditions individually"),    
                                         tableOutput("dermind"),
                                         br(),
                                         p("Now examining the CISD conditions individually, we see that **psoriasis confers a 1.2-times higher odds of having a MACE during hospitalization, and the other CISD conditions do no have a strong association, and atopic dermatitis confers a 1.3-times higher odds**, holding all other covariates constant. This is an interesting association, and so now we will analyze the CISD population in more detail."),
                                         br()
                                )
                            ),
                            tabPanel(
                                titlePanel(h4("LR MACE",style='padding-left: 15px')),
                                tabPanel(title="Pls",
                                         h3("Logistic regression for the primary outcome"),    
                                         tableOutput(outputId = "dermout"),
                                         br(),
                                         p("Now, specifically analyzing the patients with CISD, we see that beyond the non-modifiable risk factors such as age, ethnicity, and male sex, **the risk factors for MACE are prior MI (OR of 1.85), CHF (OR 6.85), COPD (OR 4.43), and interestingly, rheumatologic conditions (OR 1.43)**, when holding the other covariates constant (including admission origin, marital status, Medicare insurance, etc.). Protective factors include prior stroke, peripheral vascular disease, diabetes, and HIV. 

This has important implications for which CISD patients should be targeted for more aggressive monitoring when admitted to the hospital."),
                                         br()
                                )
                            )
                            
                            )
                        ), 
              tabPanel("Machine Learning",
                                 tabname="ml",
                                 icon=icon("calculator"),
                                 headerPanel("Machine Learning Models"),  
                       tabsetPanel(
                           tabPanel(
                               titlePanel(h4("Knn = 4",style='padding-left: 15px')),
                               tabPanel(title="Plos",
                                        h3("Knn = 4 model"),   
                                        verbatimTextOutput("confusion4"),
                                        br(),
                                        p("The accuracy of predicting in-hospital MACE is 0.89 (95% CI, 0.87, 0.90), which is the overall proportion that is predicted correctly. This can be considered good.

Based on the confusion matrix, the sensitivity (true positives divided by the sum of true positives and false negatives) is 167/(167+175) = 0.49. Thus, the false-negative rate (type-II error) is 1-sensitivity = 0.51. This means, that 51% of those who actually die are wrongly predicted to not die.

Based on the confusion matrix, the specificity (true negatives divided by the sum of true negatives and false positives) is 1799/(1799+78) = 0.96. Thus, the false-positive rate (type-I error) is 1-specificity = 0.04. This means that 4% of those who do not die are wrongly predicted to die."),
                               br()
                                        )
                           ),
                           tabPanel(
                               titlePanel(h4("Knn = 12",style='padding-left: 15px')),
                               tabPanel(title="Plos",
                                        h3("Knn = 12 model"),   
                                        verbatimTextOutput("confusion12"),
                                        br(),
                                        p("The accuracy of predicting in-hospital MACE is 0.87 (95% CI, 0.86, 0.88), which is the overall proportion that is predicted correctly. This can be considered good.

Based on the confusion matrix, the sensitivity (true positives divided by the sum of true positives and false negatives) is 108/(108+234) = 0.32. Thus, the false-negative rate (type-II error) is 1-sensitivity = 0.68. This means, that 68% of those who actually die are wrongly predicted to not die.

Based on the confusion matrix, the specificity (true negatives divided by the sum of true negatives and false positives) is 1823/(1823+54) = 0.97. Thus, the false-positive rate (type-I error) is 1-specificity = 0.03. This means that 3% of those who do not die are wrongly predicted to die."),
                                        br()
                                        
                               )
                           ),
                           tabPanel(
                               titlePanel(h4("LDA",style='padding-left: 15px')),
                               tabPanel(title="Plos",
                                        h3("LDA Model"),   
                                        verbatimTextOutput("lda"),
                                        br(),
                                        p("The accuracy of predicting in-hospital MACE is 0.85 (95% CI, 0.84, 0.86), which is the overall proportion that is predicted correctly. This can be considered good.

Based on the confusion matrix, the sensitivity (true positives divided by the sum of true positives and false negatives) is 148/(148+194) = 0.43. Thus, the false-negative rate (type-II error) is 1-sensitivity = 0.57. This means, that 57% of those who actually die are wrongly predicted to not die.

Based on the confusion matrix, the specificity (true negatives divided by the sum of true negatives and false positives) is 1737/(1737+140) = 0.93. Thus, the false-positive rate (type-I error) is 1-specificity = 0.07. This means that 7% of those who do not die are wrongly predicted to die."),
                                        br()
                               )
                           ),
                           tabPanel(
                               titlePanel(h4("RF",style='padding-left: 15px')),
                               tabPanel(title="Plos",
                                        h3("Random forest model"),   
                                        verbatimTextOutput("randomforest"),
                                        br(),
                                        p("The accuracy of predicting in-hospital MACE is 0.88 (95% CI, 0.87, 0.90), which is the overall proportion that is predicted correctly. This can be considered good.

Based on the confusion matrix, the sensitivity (true positives divided by the sum of true positives and false negatives) is 178/(178+164) = 0.52. Thus, the false-negative rate (type-II error) is 1-sensitivity = 0.48. This means, that 48% of those who actually die are wrongly predicted to not die.

Based on the confusion matrix, the specificity (true negatives divided by the sum of true negatives and false positives) is 1782/(1782+95) = 0.95. Thus, the false-positive rate (type-I error) is 1-specificity = 0.05. This means that 5% of those who do not die are wrongly predicted to die."),
                                        br()
                               )
                           ),
                           tabPanel(
                               titlePanel(h4("ROC",style='padding-left: 15px')),
                               tabPanel(title="Plos",
                                        h3("ROC curves"),   
                                        plotOutput("ROC"),
                                        br()
                                                                    )
                           ),
                           tabPanel(
                               titlePanel(h4("AUC",style='padding-left: 15px')),
                               tabPanel(title="Plos",
                                        h3("Area under the curve"),   
                                        verbatimTextOutput("auc"),
                                        br(),
                                        p("For the prediction of in-hospital MACE in patients with CISD, the best model in terms of overall accuracy and the receiver operating characteristic curve (AUROC) was the Random Forest (accuracy = 0.88 and AUC=0.89) compared to kNN with various k, and LDA. The AUC can be seen as a trade-off measure between sensitivity and the false positive rate (1-specificity). Hence, the Random Forest appears to be overall the best compromise between a relatively high sensitivity and at the same time a relatively small false positive rate, compared to the other to models. In conclusion, the Random Forest performed the best."),
                                        br()
                               )
                           )
                           
                       )

               

    )
)                
)

        
server <- function(input, output) {
    
    output$ethnicityplot <- renderPlot({
        
        unique(dermptscomplete$ethnicity)
        race <- c("WHITE", "BLACK/AFRICAN AMERICAN", "HISPANIC/LATINO", "ASIAN", "AMERICAN INDIAN/ALASKA NATIVE", "OTHER")
        class(dermptscomplete$ethnicity)
        
        dermptscomplete %>% 
            filter((ethnicity %in% race)) %>%
            dplyr::count(ethnicity = factor(ethnicity)) %>%
            mutate(pct = (prop.table(n))) %>%
            ggplot(aes(x = ethnicity, 
                       y = pct, 
                       label =  scales::percent(pct))) +
            geom_col(position = "dodge", 
                     stat = "identity", 
                     fill = "red4",
                     color = "black") +
            geom_text(position = position_dodge(width = 0.9), #move to the center of bars
                      vjust = -0.1,
                      size = 2.5) +
            scale_y_continuous(labels = scales::percent) +
            scale_x_discrete(limits = c("AMERICAN INDIAN/ALASKA NATIVE", "ASIAN", 
                                        "OTHER", "HISPANIC/LATINO", "BLACK/AFRICAN AMERICAN", "WHITE"),
                             labels = function(ethnicity) str_wrap(ethnicity, width = 
                                                                       10)) +
            xlab("Race/Ethnicity") +
            ylab("Percentages") +
            ggtitle("Race/Ethnicity of Study Population") +
            theme_light() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    output$agegenderplot = renderPlot({
        my_colors <- c("indianred1", "red4")
        
        dermptscomplete %>% 
            ggplot(aes(gender, anchor_age, fill = gender)) +
            geom_violin(width = 1, trim = FALSE) +
            geom_boxplot(width = 0.1, color = "black", alpha = 0.2) +
            scale_fill_manual(values = my_colors) +
            xlab("Gender") +
            ylab("Age") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none",
                  plot.title = element_text(size = 11)) +
            theme_minimal() +
            ggtitle("Age Distribution by Gender")
    })
    output$adminsplot = renderPlot({
        class(dermptscomplete$admission_type)
        dermptscomplete$admission_type <- as.factor(dermptscomplete$admission_type)
        unique(dermptscomplete$admission_type)
        
        dermptscomplete %>% 
            ggplot(aes(admission_type, fill = insurance)) +
            geom_bar(position = position_dodge(), color = "black") +
            scale_x_discrete(limits = c("AMBULATORY OBSERVATION", "ELECTIVE", "DIRECT OBSERVATION", "URGENT", "SURGICAL SAME DAY ADMISSION", "DIRECT EMER.", "OBSERVATION ADMIT", "EU OBSERVATION", "EW EMER."),
                             labels = function(admission_type) str_wrap(admission_type, width = 
                                                                            15)) +
            scale_fill_manual(values = c("rosybrown1", "indianred1", "red3")) +
            xlab("Admission Type") +
            ylab("Count") +
            ggtitle("Admission Classification and Insurance Type \n in the Study Population") +
            theme_light() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    output$maceethnicityplot = renderPlot({
        dermptscomplete %>%
            pivot_longer(15:41, names_to = "outcomes", values_to = "presence") %>%
            mutate(outcomes = recode(outcomes, mace_MI = "MI", mace_cardiacarrest = "Cardiac Arrest", mace_ahf = "Acute Heart Failure")) %>%
            filter(outcomes %in% c("MI", "Cardiac Arrest", "Acute Heart Failure"))%>%
            filter(presence == 1) %>%
            filter(ethnicity %in% c("AMERICAN INDIAN/ALASKA NATIVE", "ASIAN", "BLACK/AFRICAN AMERICAN", "HISPANIC/LATINO", "OTHER", "WHITE")) %>%
            ggplot(aes(x = outcomes,
                       y = anchor_age)) +
            geom_boxplot(aes(fill = ethnicity), width = 0.75) +
            scale_x_discrete(limits = c("Cardiac Arrest", "MI", 
                                        "Acute Heart Failure")) +
            xlab("Mace Outcomes") +
            ylab("Age") +
            ggtitle("Mace Outcomes by Age and Ethnicity \n in the Study Population") +
            theme_light()
    })
    output$macegenderplot <- renderPlotly({
        dermptscomplete %>%
            pivot_longer(15:41, names_to = "outcomes", values_to = "presence") %>%
            mutate(outcomes = recode(outcomes, mace_MI = "MI", mace_cardiacarrest = "Cardiac Arrest", mace_ahf = "Acute Heart Failure")) %>%
            filter(outcomes %in% c("MI", "Cardiac Arrest", "Acute Heart Failure"))%>%
            filter(presence == 1) %>%
            group_by(outcomes, gender) %>%
            ggplot(aes(x = outcomes)) +
            geom_bar(aes(fill = gender), position = position_dodge(), color = "black") +
            scale_x_discrete(limits = c("Cardiac Arrest", "MI", 
                                        "Acute Heart Failure")) +
            scale_fill_manual(values = c("indianred1", "red3")) +
            xlab("Mace Outcomes") +
            ylab("Count") +
            ggtitle("Mace Outcomes by Gender In in the Study Population") +
            theme_light()
    })
    output$macadminsplot <- renderPlotly({
        admsubjcomplete %>%
            pivot_longer(15:41, names_to = "outcomes", values_to = "presence") %>%
            mutate(outcomes = recode(outcomes, mace_MI = "MI", mace_cardiacarrest = "Cardiac Arrest", mace_ahf = "Acute Heart Failure")) %>%
            filter(outcomes %in% c("MI", "Cardiac Arrest", "Acute Heart Failure"))%>%
            filter(presence == 1) %>%
            ggplot(aes(x = outcomes,
                       y = anchor_age)) +
            geom_boxplot(aes(fill = ethnicity), width = 0.75) +
            scale_x_discrete(limits = c("Cardiac Arrest", "MI", 
                                        "Acute Heart Failure")) +
            xlab("Mace Outcomes") +
            ylab("Age") +
            ggtitle("Mace Outcomes by Age and Ethnicity \n in the Study Population") +
            theme_light()
    })    
    
    output$dermwhole <- function(){
        knitr::kable(tab) %>%
            kable_styling("striped", full_width = F) 
    }   
    output$dermind <- function(){
        knitr::kable(tab2) %>%
            kable_styling("striped", full_width = F) 
    }
    output$dermout <- function(){
        knitr::kable(tab3) %>%
            kable_styling("striped", full_width = F)
    }    
    
    output$confusion4 <- renderPrint({
        data <- dermptscomplete %>% dplyr::select(subject_id,mace, gender, age_cat, ethnicity, mi, chf, hiv, pvd, cva, copd, rheum,  livmodsev, malig, miliv, pleg, dm, numadm) %>% 
            mutate(gender = as.factor(gender), 
                   age_cat = as.factor(age_cat), 
                   ethnicity = as.factor(ethnicity), 
                   mi = as.factor(mi), 
                   chf = as.factor(chf),
                   hiv = as.factor(hiv), 
                   pvd = as.factor(pvd), 
                   cva = as.factor(cva), 
                   copd = as.factor(copd), 
                   rheum = as.factor(rheum), 
                   livmodsev = as.factor(livmodsev),
                   malig = as.factor(malig),
                   miliv = as.factor(miliv),
                   pleg = as.factor(pleg), 
                   dm = as.factor(dm),
                   numadm = as.numeric(numadm),
                   mace = as.factor(mace)) 

        # Use complete cases 
        data <- data[complete.cases(data), ]
        
        set.seed(1)
        
        x <- stratified(data, "mace", 0.7, keep.rownames = TRUE)
        train_set <- x %>% dplyr::select(-rn)
        train_index <- as.numeric(x$rn)
        test_set <- data[-train_index,]
        
        ###### KNN with k=4
        knn_fit <- knn3(mace~.,data = dplyr::select(train_set, mace, gender, age_cat, ethnicity, mi, chf, hiv, pvd, cva, copd, rheum, livmodsev, malig, miliv, pleg, dm, numadm), k = 4)
        f_hat_4 <- predict(knn_fit, newdata = test_set)[,2]
        tab_4 <- table(pred=round(f_hat_4), truth=test_set$mace)
        confusionMatrix(tab_4, positive = "1")
    })
    
    output$confusion12 <- renderPrint({
        data <- dermptscomplete %>% dplyr::select(subject_id,mace, gender, age_cat, ethnicity, mi, chf, hiv, pvd, cva, copd, rheum,  livmodsev, malig, miliv, pleg, dm, numadm) %>% 
            mutate(gender = as.factor(gender), 
                   age_cat = as.factor(age_cat), 
                   ethnicity = as.factor(ethnicity), 
                   mi = as.factor(mi), 
                   chf = as.factor(chf),
                   hiv = as.factor(hiv), 
                   pvd = as.factor(pvd), 
                   cva = as.factor(cva), 
                   copd = as.factor(copd), 
                   rheum = as.factor(rheum), 
                   livmodsev = as.factor(livmodsev),
                   malig = as.factor(malig),
                   miliv = as.factor(miliv),
                   pleg = as.factor(pleg), 
                   dm = as.factor(dm),
                   numadm = as.numeric(numadm),
                   mace = as.factor(mace)) 
        
        # Use complete cases 
        data <- data[complete.cases(data), ]
        
        set.seed(1)
        
        x <- stratified(data, "mace", 0.7, keep.rownames = TRUE)
        train_set <- x %>% dplyr::select(-rn)
        train_index <- as.numeric(x$rn)
        test_set <- data[-train_index,]
        
        knn_fit <- knn3(mace~.,data = dplyr::select(train_set, mace, gender, age_cat, ethnicity, mi, chf, hiv, pvd, cva, copd, rheum, livmodsev, malig, miliv, pleg, dm, numadm), k =12)
        f_hat_12 <- predict(knn_fit, newdata = test_set)[,2]
        tab_12 <- table(pred=round(f_hat_12), truth=test_set$mace)
        confusionMatrix(tab_12, positive = "1")
    })
    
    output$lda <- renderPrint({
        data <- dermptscomplete %>% dplyr::select(subject_id,mace, gender, age_cat, ethnicity, mi, chf, hiv, pvd, cva, copd, rheum,  livmodsev, malig, miliv, pleg, dm, numadm) %>% 
            mutate(gender = as.factor(gender), 
                   age_cat = as.factor(age_cat), 
                   ethnicity = as.factor(ethnicity), 
                   mi = as.factor(mi), 
                   chf = as.factor(chf),
                   hiv = as.factor(hiv), 
                   pvd = as.factor(pvd), 
                   cva = as.factor(cva), 
                   copd = as.factor(copd), 
                   rheum = as.factor(rheum), 
                   livmodsev = as.factor(livmodsev),
                   malig = as.factor(malig),
                   miliv = as.factor(miliv),
                   pleg = as.factor(pleg), 
                   dm = as.factor(dm),
                   numadm = as.numeric(numadm),
                   mace = as.factor(mace)) 
        
        # Use complete cases 
        data <- data[complete.cases(data), ]
        
        set.seed(1)
        
        x <- stratified(data, "mace", 0.7, keep.rownames = TRUE)
        train_set <- x %>% dplyr::select(-rn)
        train_index <- as.numeric(x$rn)
        test_set <- data[-train_index,]
        
        set.seed(1)
        
        lda_fit <- lda(mace ~ gender + age_cat + ethnicity + mi + chf + hiv + pvd + cva + copd + rheum + livmodsev + malig + miliv + pleg + dm + numadm, data = train_set)
        
        p_hat_lda <- predict(lda_fit, newdata = test_set)$posterior[,2]
        
        y_hat_lda <- factor(ifelse(p_hat_lda > 0.5, 1, 0))
        
        confusionMatrix(data = as.factor(y_hat_lda), reference = test_set$mace, positive = "1")
    })
    output$randomforest <- renderPrint({
        data <- dermptscomplete %>% dplyr::select(subject_id,mace, gender, age_cat, ethnicity, mi, chf, hiv, pvd, cva, copd, rheum,  livmodsev, malig, miliv, pleg, dm, numadm) %>% 
            mutate(gender = as.factor(gender), 
                   age_cat = as.factor(age_cat), 
                   ethnicity = as.factor(ethnicity), 
                   mi = as.factor(mi), 
                   chf = as.factor(chf),
                   hiv = as.factor(hiv), 
                   pvd = as.factor(pvd), 
                   cva = as.factor(cva), 
                   copd = as.factor(copd), 
                   rheum = as.factor(rheum), 
                   livmodsev = as.factor(livmodsev),
                   malig = as.factor(malig),
                   miliv = as.factor(miliv),
                   pleg = as.factor(pleg), 
                   dm = as.factor(dm),
                   numadm = as.numeric(numadm),
                   mace = as.factor(mace)) 
        
        # Use complete cases 
        data <- data[complete.cases(data), ]
        
        set.seed(1)
        
        x <- stratified(data, "mace", 0.7, keep.rownames = TRUE)
        train_set <- x %>% dplyr::select(-rn)
        train_index <- as.numeric(x$rn)
        test_set <- data[-train_index,]
        
        RNGkind(sample.kind = "Rounding")
        set.seed(1)
        
        forest_fit <- randomForest(mace ~ gender + age_cat + ethnicity + mi + chf + hiv + pvd + cva + copd + rheum + livmodsev + malig + miliv + pleg + dm + numadm, data = train_set)
        
        p_hat_forest <- predict(forest_fit, newdata = test_set, type = "prob")[, 2]
        
        y_hat_forest <- factor(ifelse(p_hat_forest > 0.5, 1, 0))
        
        confusionMatrix(data = as.factor(y_hat_forest), reference = test_set$mace, positive = "1")
    })
    output$ROC = renderPlot({
        knn_4_roc <- roc(test_set$mace, f_hat_4)
        knn_8_roc <- roc(test_set$mace, f_hat_8)
        knn_12_roc <- roc(test_set$mace, f_hat_12)
        roc_lda <- roc(test_set$mace, p_hat_lda)
        roc_forest <- roc(test_set$mace, p_hat_forest)
        
        ggroc(list("kNN, k=4" = knn_4_roc, "kNN, k=8" = knn_8_roc, "kNN, k=12" = knn_12_roc, "LDA" = roc_lda, "Random Forest" = roc_forest), legacy.axes = F)+
            geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "black", linetype = "dashed") +
            ylab("Sensitivity") +
            xlab("Specificity") +
            labs(color = "Models") + 
            ggtitle("ROC Curves for 5 Models") +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5))
    })
    output$auc <- renderPrint({
   AUCS
    })
    } 

shinyApp(ui = ui, server = server)
