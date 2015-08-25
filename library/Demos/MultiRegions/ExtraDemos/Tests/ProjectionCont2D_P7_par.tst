<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>C0 Projection 2D with P=7, parallel</description>
    <executable>ProjectCont2D</executable>
    <parameters>-I GlobalSysSoln=IterativeStaticCond -P IterativeSolverTolerance=1e-12 ProjectionCont2D_P7.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">ProjectionCont2D_P7.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-6">5.79152e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-5">6.52922e-06</value>
        </metric>
    </metrics>
</test>


