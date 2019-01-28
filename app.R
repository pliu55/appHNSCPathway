#!/usr/bin/env Rscript

#
#  pliu 20190125
#
#  make a Shiny App
#
#  https://www.shinyapps.io/admin/#/dashboard
#

library(data.table)
library(shiny)
library(visNetwork)
suppressMessages(library(plotly))

library(RColorBrewer)
library(scales)

readDtList <- function(){
    fnode = './node.tsv'
    fedge = './edge.tsv'
    fmed  = './median.tsv'

    nodedt = fread(fnode, header=TRUE, sep="\t")
    edgedt = fread(fedge, header=TRUE, sep="\t")
    meddt  = fread(fmed,  header=TRUE, sep="\t")

    dtlist = list( med   = meddt, 
                   node  = nodedt,
                   edge  = edgedt )

    return(dtlist)
}


dtlist <- readDtList()


ui <- fluidPage(
    titlePanel( 'HNSC NCI pathways' ),
    hr(),


    fluidRow(
        column(width = 3, 
               plotlyOutput('pointPlot', width='260px', height='260px'),
               br(),
               selectizeInput(  inputId  = 'pathway', 
                                label    = 'select/enter a pathway name below:',
                                multiple = FALSE,
                                choices  = NULL ),
               br(),
               br(),
               htmlOutput( outputId = 'disclaimer' )
        ),

        column(width  = 4, 
               offset = 0.1,
               visNetworkOutput( outputId = 'ntw_pos',
                                 height   = '400px'  ),
             # visNetworkOutput( outputId = 'ntw_shape',
             #                   width    = '300px',
             #                   height   = '80px' )
             # textOutput( outputId = 'legend')
               tableOutput('legend_shape')
        ),
        
        column(width  = 4,
               offset = 0.1,
               visNetworkOutput( outputId = 'ntw_neg',
                                 height   = '400px'  ),
             # visNetworkOutput( outputId = 'ntw_color',
             #                   width    = '200px',
             #                   height   = '80px' )
             # textOutput('text_neg')
             # tableOutput('neg_nnode')
               tableOutput('legend_color')
        )
    )
)


server <- function(input, output, session) {
    meddt = dtlist$med 
    pos_nodedt = reactive(dtlist$node[ (hpv == 'pos') & (pth == input$pathway)])
    neg_nodedt = reactive(dtlist$node[ (hpv == 'neg') & (pth == input$pathway)])
    edgedt = reactive(dtlist$edge[ pth == input$pathway ])

    updateSelectizeInput( session, 
                          inputId  = 'pathway', 
                          choices  = sort(meddt$name), 
                          selected = 'RhoA signaling pathway',
                          server   = TRUE )

    seldt = reactive(meddt[ name %in% input$pathway])

    output$pointPlot = renderPlotly(makePlotly(meddt, seldt, input))

    output$ntw_pos = renderVisNetwork(plotNtw(pos_nodedt, edgedt, input))
    output$ntw_neg = renderVisNetwork(plotNtw(neg_nodedt, edgedt, input))

    output$legend_color = renderTable({
        data.table( 
            color = c( '<font color="#FB6A4A"><strong>red</strong></font>',
                       '<font color="#6BAED6"><strong>blue</strong></font>',
                       '<font color="#D9D9D9"><strong>grey</strong></font>' ),
            `correlated with overall survival (p < 0.05)` = c(
                'both pathway level and geomic data',
                'only pathway level',
                'none' ) )
    }, align='c', sanitize.text.function = function(x) x)

    output$legend_shape = renderTable({
        data.table( shape = c( '&#9733;', '&#9650;', '&#11044;' ),
            `patients with perturbed pathway level` = c(
                '> 50%', '&#8804; 50%', 'none' ))
    }, align='c', sanitize.text.function = function(x) x)

    output$disclaimer = renderText({
        paste0('<font size="2px">',
               'This app is an ongoing work supported in part by <a href="https://hn-spore.wisc.edu/cep/awarded-pilots/#gitter">a Wisconsin Head and Neck Cancer SPORE Career Enhancement Award</a> to Anthony Gitter in collaboration with Paul Ahlquist, David Page, Irene Ong, and Peng Liu.',
               '<br><br>',
               'For question or feature request, please create an issue at this app\'s <a href="https://github.com/pliu55/appHNSCPathway/issues">GitHub page</a>.',
               '</font>')
    })
}


plotNtw <- function(nodedt, edgedt, input) {
    tlist = list( text  = '', 
                  style = paste0('font-family:Arial;font-weight:normal;',
                                 'font-size:17px;text-align:center') )
    if ( nodedt()$hpv[1] %in% c( 'pos') ){ 
        tlist$text = 'HPV(+)'
    } else if ( nodedt()$hpv[1] %in% c( 'neg') ) {
        tlist$text = 'HPV(-)'
    }

    visNetwork(nodes = nodedt(), edges = edgedt(), main=tlist) %>%
        visOptions(highlightNearest = list(enabled=TRUE, degree=0)) %>%
        visLayout(randomSeed=12345)
}


makePlotly <- function(meddt, seldt, input) {
    p = ggplot(meddt, aes(x=`HPV(+)`, y=`HPV(-)`, text=name)) + 
        geom_abline(slope=1, intercept=0, linetype='dashed', size=0.3) +
        geom_point(size=2, shape=1) + 
        geom_point(data=seldt(), aes(x=`HPV(+)`, y=`HPV(-)`),
                   color='green2', size=2, shape=1) +
        theme_bw() + 
        theme( plot.title = element_text(hjust=0.5), 
               axis.text  = element_text(size=9) ) +
        ggtitle( 'pathway perturbation' ) +
        xlab( 'HPV(+) patients' ) +
        ylab( 'HPV(-) patients' )
        
    ply = ggplotly(p, tooltip=c('text'))
    ply$elementId = NULL
    layout( ply,
            height   = input$plotHeight, 
            autosize = TRUE)
}


shinyApp(ui, server)

