<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>HDG Helmholtz 2D CG with P=7</description>
    <executable>HDGHelmholtz2D</executable>
    <parameters>Helmholtz2D_P7_AllBCs.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Helmholtz2D_P7_AllBCs.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">0.011251</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.0355625</value>
        </metric>
    </metrics>
</test>
