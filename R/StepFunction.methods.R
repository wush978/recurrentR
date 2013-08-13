evalqOnLoad({
	
setMethod("^",
	signature(e1 = StepFunction),
	function (e1, e2) 
	{
		new(StepFunction, e1$x, e1$y^e2)
	}
)
	
	
})