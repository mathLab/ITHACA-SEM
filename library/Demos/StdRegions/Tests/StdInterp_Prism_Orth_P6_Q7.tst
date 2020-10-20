<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpBasis Prism Orthonormal basis P=6 Q=7</description>
    <executable>StdInterp</executable>
    <parameters>-s prism -b Ortho_A Ortho_A Ortho_B -o 6 6 6 -p 7 7 7 -P GaussGaussLegendre GaussGaussLegendre GaussGaussLegendre</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.89694e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.77636e-14</value>
        </metric>
    </metrics>
</test>

