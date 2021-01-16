<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterp Seg Modified_A basis P=7 Q=8</description>
    <executable>StdInterp</executable>
    <parameters>-s segment -b Modified_A -o 7 -p 8 -P GaussGaussLegendre</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">2.81216e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.33227e-15</value>
        </metric>
    </metrics>
</test>
