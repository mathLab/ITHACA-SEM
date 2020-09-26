<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterp Hex Modified basis P=6 Q=7</description>
    <executable>StdInterp</executable>
    <parameters>-s hexahedron -b Modified_A Modified_A Modified_A -o 6 6 6 -p 7
7 7 -P GaussGaussLegendre GaussGaussLegendre GaussGaussLegendre</parameters>
    <metrics>
        <metric type="Linf" id="1">
            <value tolerance="1e-12">0</value>
        </metric>
        <metric type="L2" id="2">
            <value tolerance="1e-12">0</value>
        </metric>
    </metrics>
</test>


