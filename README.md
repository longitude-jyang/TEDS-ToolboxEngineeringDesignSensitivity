# Toolbox for Engineering Design Sensitivity (TEDS)
Research relating to Digital Twins (DTs) has mostly focused on manufacturing, operation and maintenance, where the corre sponding physical entities exist. However, at the design stage, the virtual part of a DT predates its physical part. As a result, the DT at the design stage is simulation based and mirrors the expected functions and events of the physical entity. With an emphasis on data integration, an overview of the application of a DT in the design process under uncertainty is given below. At each step of the design process, the quantities of interest (QoIs) are estimated from the DT. Decision makers then extract useful information form QoIs using design metrics and decide either to accept the design, to improve it, or to obtain more data to better inform the decision making process.

![Alt Text](/docs/singleslide_process.gif)

As the design requirements are evolving during the design process, we have divided the design process into two stages, the early conceptual design stage and the detailed design stage. At detailed stages of the design process, we normally have a specified design target or key performance indicator (KPI). 

The need is to find the probability that the design will meet the target, $P_\text{f}$, and also the sensitivity of this probability to input data so that design decisions can be made accordingly. In this case, the design metrics are denoted as KPI-based and are most suitable for detailed design stages. 

On the other hand, at early design stages, a specific design requirement is not normally available. This could be due to either too much uncertainty to fix the requirements, or a preference to waive specific requirements so that a much wider design space can be explored. As there are no specific design KPIs, the design metrics that are most suitable at the conceptual design stage are called KPI-free metrics. The KPI-free metrics are the design entropy H and its sensitivity to perturbations of design parameters, which are independent of specific design targets (KPIs).

An introduction of TEDS and the project background is given here [Introduction of TEDS](/docs/TEDS_ToolboxEngineeringDesignSensitivity_Git.pdf) 

![TEDS for design process](/docs/TEDS_1.png) 
