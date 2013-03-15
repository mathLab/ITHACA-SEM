<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=2</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Channel_Flow_3modes_Parallel.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Channel_Flow_3modes_Parallel.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9">5.09263e-09</value>
            <value variable="v" tolerance="1e-9">1.58852e-09</value>
            <value variable="p" tolerance="1e-7">1.24401e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">4.43727e-08</value>
            <value variable="v" tolerance="1e-9">3.64497e-09</value>
            <value variable="p" tolerance="1e-7">2.71649e-07</value>
        </metric>
    </metrics>
</test>


